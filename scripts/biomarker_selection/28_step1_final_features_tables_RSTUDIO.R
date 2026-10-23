#======PARAMETROS===========

SELECTION_RUN <- path.expand(
  "~/Documents/PhD/metatest/results_reproducible_20261010/comparacion_modelamientos/upstream_comun/downstream/seleccion_mejor_escenario_train/20261011T233118Z_pid679971"
)

MIN_MODEL_VOTES <- 2L
VIP_THRESHOLD <- 1
OPEN_VIEWS <- TRUE

FAMILIES <- c("P1", "P2", "P3")
CLASS_LEVELS <- c("NP", "FD", "SO")


DIRECTION_LEVELS <- c("FD", "NP", "SO")

VIEW_LABELS <- c(
  tx = "TX",
  pr = "PR",
  me = "ME",
  cl = "CL"
)

#======PARAMETROS GRAFICOS GLOBALES===========

#----- BARRAS -----
BAR_BASE_SIZE <- 20
BAR_TEXT_SIZE <- 26
BAR_TITLE_SIZE <- 26
BAR_STRIP_TEXT_SIZE <- 21
BAR_LEGEND_TEXT_SIZE <- 17
BAR_LEGEND_TITLE_SIZE <- 18

#----- HEATMAPS -----
HM_BASE_SIZE <- 15
HM_TITLE_SIZE <- 21
HM_AXIS_X_SIZE <- 14
HM_AXIS_Y_SIZE <- 14
HM_CELL_TEXT_SIZE <- 5
HM_STRIP_TEXT_SIZE <- 14
HM_LEGEND_TEXT_SIZE <- 15
HM_LEGEND_TITLE_SIZE <- 17
HM_X_ANGLE <- 50

#----- RED ACOPLADA (coupled network) -----
NET_BASE_SIZE <- 14
NET_FEATURE_LABEL_SIZE <- 5.5
NET_PATHWAY_LABEL_SIZE <- 5.2
NET_FEATURE_NODE_SIZE_RANGE <- c(5.0, 13.0)
NET_PATHWAY_NODE_SIZE <- 6.2
NET_LEGEND_TEXT_SIZE <- 12
NET_LEGEND_TITLE_SIZE <- 13
NET_PLOT_MARGIN <- c(1.6, 3.0, 2.5, 3.0)
NET_COORD_XLIM <- c(-11.5, 11.5)
NET_COORD_YLIM <- c(-7.5, 7.5)
NET_PATH_LABEL_WRAP <- 18
NET_GGREPEL_MAX_TIME <- 6

# Etiquetas explicitas para que FM / PM se entiendan
NETWORK_FM_LABELS <- c(
  FM1 = "FM1 = feature module 1",
  FM2 = "FM2 = feature module 2",
  FM3 = "FM3 = feature module 3",
  FM4 = "FM4 = feature module 4",
  FM5 = "FM5 = feature module 5"
)

NETWORK_PM_LABELS <- c(
  PM1 = "PM1 = pathway module 1",
  PM2 = "PM2 = pathway module 2",
  PM3 = "PM3 = pathway module 3",
  PM4 = "PM4 = pathway module 4",
  PM5 = "PM5 = pathway module 5"
)

#----- RED TRIPARTITA -----
TRIP_BASE_SIZE <- 13
TRIP_FEATURE_LABEL_SIZE <- 5.6
TRIP_PROGRAM_LABEL_SIZE <- 5.8
TRIP_PATHWAY_LABEL_SIZE <- 5.4
TRIP_FEATURE_SIZE_RANGE <- c(3.2, 10.5)
TRIP_FP_EDGE_RANGE <- c(0.6, 2.4)
TRIP_PP_EDGE_RANGE <- c(1.0, 4.2)
TRIP_PROGRAM_NODE_SIZE <- 7.5
TRIP_PATHWAY_NODE_SIZE <- 4.5
TRIP_LEGEND_TEXT_SIZE <- 11
TRIP_LEGEND_TITLE_SIZE <- 12
TRIP_PLOT_MARGIN <- c(1, 8, 1, 8)

#----- MODULOS / CLUSTERS EN REDES -----
NETWORK_K_FEATURE_MODULES <- 4L
NETWORK_K_PATHWAY_MODULES <- 4L
#======paleta de aristas por programa===========

EDGE_PROGRAM_COLORS <- c(
  P01 = "#1b9e77",
  P02 = "#d95f02",
  P03 = "#7570b3"
)
#======FUNCIONES===========

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
  x <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  missing_cols <- setdiff(required, names(x))
  
  assert(
    length(missing_cols) == 0L,
    paste(
      "Faltan columnas en",
      basename(path),
      ":",
      paste(missing_cols, collapse = ", ")
    )
  )
  
  x
}

find_latest_bio_run <- function(selection_run) {
  #======resolver BIOLOGICAL_PROGRAMS===========
  root <- file.path(
    selection_run,
    "BIOLOGICAL_PROGRAMS"
  )
  
  assert(
    dir.exists(root),
    paste("No existe:", root)
  )
  
  latest_file <- file.path(
    root,
    "LATEST_RUN.txt"
  )
  
  if (file.exists(latest_file)) {
    p <- trimws(
      readLines(
        latest_file,
        warn = FALSE
      )[1L]
    )
    
    if (nzchar(p) && dir.exists(p)) {
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
    "No se encontro ninguna corrida en BIOLOGICAL_PROGRAMS."
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

find_latest_model_run <- function(
    bio_run,
    subdir,
    required_rel
) {
  #======resolver corrida de modelo===========
  root <- file.path(
    bio_run,
    subdir
  )
  
  assert(
    dir.exists(root),
    paste("No existe:", root)
  )
  
  latest_file <- file.path(
    root,
    "LATEST_RUN.txt"
  )
  
  if (file.exists(latest_file)) {
    p <- trimws(
      readLines(
        latest_file,
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
        required_rel
      )
    )
  ]
  
  assert(
    length(dirs) > 0L,
    paste(
      "No se encontro corrida completa en",
      root
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


find_crossomic_input_audit <- function(
    selection_run,
    family
) {
  #======localizar auditoria crossomic con ruta del modelo MOFA===========
  root <- file.path(
    selection_run,
    family,
    "BIDIRECTIONAL_FRAMEWORK"
  )
  
  candidates <- list.files(
    root,
    pattern = "^00_INPUT_AUDIT\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  candidates <- candidates[
    grepl(
      "/CROSSOMIC_BRIDGE/.*/tables/00_INPUT_AUDIT\\.csv$",
      candidates
    )
  ]
  
  assert(
    length(candidates) > 0L,
    paste(
      "No se encontro 00_INPUT_AUDIT.csv del CROSSOMIC_BRIDGE para",
      family
    )
  )
  
  info <- file.info(
    candidates
  )
  
  candidates[
    order(
      info$mtime,
      basename(candidates),
      decreasing = TRUE
    )
  ][1L]
}

map_mofa_view <- function(x) {
  #======traducir nombres de vistas MOFA a view_code===========
  z <- tolower(
    iconv(
      as.character(x),
      from = "",
      to = "ASCII//TRANSLIT"
    )
  )
  
  out <- rep(
    NA_character_,
    length(z)
  )
  
  out[grepl("transcript", z)] <- "tx"
  out[grepl("prote", z)] <- "pr"
  out[grepl("metabol", z)] <- "me"
  out[grepl("clinic", z)] <- "cl"
  
  assert(
    all(!is.na(out)),
    paste(
      "No se pudieron mapear estas vistas MOFA:",
      paste(
        unique(x[is.na(out)]),
        collapse = ", "
      )
    )
  )
  
  out
}

extract_mofa_importance <- function(
    model_path,
    family
) {
  #======extraer importancia nativa MOFA del modelo final===========
  assert(
    file.exists(model_path),
    paste(
      "No existe el modelo MOFA:",
      model_path
    )
  )
  
  assert(
    requireNamespace(
      "MOFA2",
      quietly = TRUE
    ),
    "Se necesita el paquete MOFA2."
  )
  
  model <- readRDS(
    model_path
  )
  
  weights <- MOFA2::get_weights(
    model,
    views = "all",
    factors = "all",
    as.data.frame = TRUE
  )
  
  assert(
    all(
      c(
        "feature",
        "factor",
        "value",
        "view"
      ) %in%
        names(weights)
    ),
    "get_weights() no devolvio las columnas esperadas."
  )
  
  weights$family <- family
  
  weights$view_code <- map_mofa_view(
    weights$view
  )
  
  weights$feature_original <- as.character(
    weights$feature
  )
  
  weights$value <- suppressWarnings(
    as.numeric(
      weights$value
    )
  )
  
  weights$loading_squared <-
    weights$value^2
  
  summary <- stats::aggregate(
    loading_squared ~
      family +
      view_code +
      feature_original,
    data = weights,
    FUN = sum,
    na.rm = TRUE
  )
  
  summary$MOFA_loading_norm <- sqrt(
    summary$loading_squared
  )
  
  summary$loading_squared <- NULL
  
  summary
}

percentile_rank <- function(x) {
  #======normalizar una metrica a percentil 0-1===========
  out <- rep(
    NA_real_,
    length(x)
  )
  
  ok <- is.finite(
    x
  )
  
  if (!any(ok)) {
    return(out)
  }
  
  values <- x[ok]
  
  if (
    length(values) == 1L ||
    length(unique(values)) == 1L
  ) {
    out[ok] <- 0.5
    return(out)
  }
  
  out[ok] <- (
    rank(
      values,
      ties.method = "average"
    ) - 1
  ) / (
    length(values) - 1
  )
  
  out
}

framework_pls_margin <- function(
    direction,
    align_np,
    align_fd,
    align_so
) {
  #======margen PLS respecto a la direccion fijada por el framework===========
  out <- rep(
    NA_real_,
    length(direction)
  )
  
  for (i in seq_along(direction)) {
    d <- direction[i]
    
    a <- c(
      NP = align_np[i],
      FD = align_fd[i],
      SO = align_so[i]
    )
    
    if (
      is.na(d) ||
      !d %in% names(a) ||
      !is.finite(a[d])
    ) {
      next
    }
    
    other <- a[
      names(a) != d
    ]
    
    other <- other[
      is.finite(other)
    ]
    
    if (!length(other)) {
      next
    }
    
    out[i] <- a[d] - max(
      other
    )
  }
  
  out
}

pattern_to_class <- function(x) {
  #======traducir patron a lado fenotipico===========
  out <- rep(
    NA_character_,
    length(x)
  )
  
  out[x == "NP_like_high"] <- "NP"
  out[x == "FD_like_high"] <- "FD"
  out[x == "SO_like_high"] <- "SO"
  
  out
}

build_pls_direction <- function(
    pls17_run,
    pls_loadings,
    pls_status
) {
  #======reconstruir direccion PLS desde resultados congelados===========
  out <- list()
  
  for (fam in FAMILIES) {
    chosen_ncomp <- suppressWarnings(
      as.integer(
        pls_status$chosen_ncomp[
          match(
            fam,
            pls_status$family
          )
        ]
      )
    )
    
    assert(
      is.finite(chosen_ncomp) &&
        chosen_ncomp >= 1L,
      paste(
        "chosen_ncomp invalido para",
        fam
      )
    )
    
    score_path <- file.path(
      pls17_run,
      "tables",
      fam,
      paste0(
        fam,
        "_08_PLSDA_SCORES.csv"
      )
    )
    
    scores <- read_csv_checked(
      score_path,
      required = c(
        "grupo",
        "split"
      )
    )
    
    score_cols <- paste0(
      "PLS",
      seq_len(chosen_ncomp)
    )
    
    assert(
      all(
        score_cols %in%
          names(scores)
      ),
      paste(
        "Faltan scores PLS en",
        fam
      )
    )
    
    train_scores <- scores[
      scores$split == "TRAIN",
      ,
      drop = FALSE
    ]
    
    class_centroids <- do.call(
      rbind,
      lapply(
        CLASS_LEVELS,
        function(cls) {
          z <- train_scores[
            train_scores$grupo == cls,
            score_cols,
            drop = FALSE
          ]
          
          data.frame(
            grupo = cls,
            t(
              colMeans(
                z,
                na.rm = TRUE
              )
            ),
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
        }
      )
    )
    
    load_fam <- pls_loadings[
      pls_loadings$family == fam &
        pls_loadings$component <= chosen_ncomp,
      ,
      drop = FALSE
    ]
    
    rows <- lapply(
      split(
        load_fam,
        load_fam$predictor_key
      ),
      function(z) {
        beta <- rep(
          0,
          chosen_ncomp
        )
        
        for (cc in seq_len(chosen_ncomp)) {
          hit <- z$loading[
            z$component == cc
          ]
          
          if (length(hit)) {
            beta[cc] <- hit[1L]
          }
        }
        
        alignment <- vapply(
          CLASS_LEVELS,
          function(cls) {
            mu <- as.numeric(
              class_centroids[
                class_centroids$grupo == cls,
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
          alignment,
          decreasing = TRUE
        )
        
        data.frame(
          family = fam,
          predictor_key = z$predictor_key[1L],
          PLS_direction = names(alignment)[ord[1L]],
          PLS_direction_margin = if (
            length(ord) >= 2L
          ) {
            unname(
              alignment[ord[1L]] -
                alignment[ord[2L]]
            )
          } else {
            NA_real_
          },
          PLS_alignment_NP = unname(
            alignment["NP"]
          ),
          PLS_alignment_FD = unname(
            alignment["FD"]
          ),
          PLS_alignment_SO = unname(
            alignment["SO"]
          ),
          stringsAsFactors = FALSE
        )
      }
    )
    
    out[[fam]] <- do.call(
      rbind,
      rows
    )
  }
  
  final <- do.call(
    rbind,
    out
  )
  
  rownames(final) <- NULL
  
  final
}

build_sharing_table <- function(all_features) {
  #======comparar P1 P2 P3 por feature===========
  ids <- unique(
    all_features$feature_id
  )
  
  rows <- lapply(
    ids,
    function(id) {
      z <- all_features[
        all_features$feature_id == id,
        ,
        drop = FALSE
      ]
      
      get_value <- function(fam, col) {
        hit <- z[
          z$family == fam,
          col
        ]
        
        if (!length(hit)) {
          return(NA)
        }
        
        hit[1L]
      }
      
      present <- FAMILIES[
        FAMILIES %in%
          z$family
      ]
      
      data.frame(
        feature_id = id,
        view_code = z$view_code[1L],
        feature_original = z$feature_original[1L],
        feature_label = z$feature_label[1L],
        
        in_P1 = "P1" %in% z$family,
        in_P2 = "P2" %in% z$family,
        in_P3 = "P3" %in% z$family,
        
        support_P1 = get_value(
          "P1",
          "n_model_votes"
        ),
        support_P2 = get_value(
          "P2",
          "n_model_votes"
        ),
        support_P3 = get_value(
          "P3",
          "n_model_votes"
        ),
        
        tier_P1 = get_value(
          "P1",
          "support_tier"
        ),
        tier_P2 = get_value(
          "P2",
          "support_tier"
        ),
        tier_P3 = get_value(
          "P3",
          "support_tier"
        ),
        
        score_P1 = get_value(
          "P1",
          "evidence_rank_score"
        ),
        score_P2 = get_value(
          "P2",
          "evidence_rank_score"
        ),
        score_P3 = get_value(
          "P3",
          "evidence_rank_score"
        ),
        
        pattern_P1 = get_value(
          "P1",
          "framework_pattern"
        ),
        pattern_P2 = get_value(
          "P2",
          "framework_pattern"
        ),
        pattern_P3 = get_value(
          "P3",
          "framework_pattern"
        ),
        
        n_questions = length(
          present
        ),
        questions = paste(
          present,
          collapse = "+"
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  
  out <- do.call(
    rbind,
    rows
  )
  
  out <- out[
    order(
      -out$n_questions,
      out$view_code,
      out$feature_label
    ),
    ,
    drop = FALSE
  ]
  
  rownames(out) <- NULL
  
  out
}

#======parametros piloto modulos feature-funcion===========
PILOT_FAMILY <- "P1"
PILOT_PATTERN <- "SO_like_high"
PILOT_N_FEATURE_MODULES <- 3L
PILOT_N_PATHWAY_MODULES <- 3L
PILOT_CLUSTER_METHOD <- "average"

#======funcion matriz simetrica desde correlaciones long===========
build_symmetric_corr_matrix <- function(
    corr_long,
    features,
    feature1_col = "feature_1",
    feature2_col = "feature_2",
    value_col = "rho"
) {
  
  features <- unique(as.character(features))
  n <- length(features)
  
  out <- matrix(
    NA_real_,
    nrow = n,
    ncol = n,
    dimnames = list(features, features)
  )
  
  diag(out) <- 1
  
  if (nrow(corr_long) > 0L) {
    for (i in seq_len(nrow(corr_long))) {
      f1 <- as.character(corr_long[[feature1_col]][i])
      f2 <- as.character(corr_long[[feature2_col]][i])
      vv <- suppressWarnings(as.numeric(corr_long[[value_col]][i]))
      
      if (is.na(vv)) next
      if (!(f1 %in% features) || !(f2 %in% features)) next
      
      out[f1, f2] <- vv
      out[f2, f1] <- vv
    }
  }
  
  out
}

#======funcion matriz feature x pathway===========
build_feature_pathway_matrix <- function(
    assoc_long,
    feature_col = "feature",
    pathway_col = "pathway",
    value_col = "association_percentile"
) {
  
  df <- assoc_long
  
  df[[feature_col]] <- as.character(df[[feature_col]])
  df[[pathway_col]] <- as.character(df[[pathway_col]])
  df[[value_col]] <- suppressWarnings(as.numeric(df[[value_col]]))
  df[[value_col]][is.na(df[[value_col]])] <- 0
  
  feats <- unique(df[[feature_col]])
  paths <- unique(df[[pathway_col]])
  
  mat <- matrix(
    0,
    nrow = length(feats),
    ncol = length(paths),
    dimnames = list(feats, paths)
  )
  
  for (i in seq_len(nrow(df))) {
    rr <- df[[feature_col]][i]
    cc <- df[[pathway_col]][i]
    vv <- df[[value_col]][i]
    mat[rr, cc] <- vv
  }
  
  mat
}

#======funcion clustering de modulos===========
build_feature_function_modules <- function(
    assoc_long,
    corr_long,
    family,
    pattern,
    n_feature_modules = 3L,
    n_pathway_modules = 3L,
    feature_col = "feature",
    pathway_col = "pathway",
    assoc_col = "association_percentile",
    corr_feature1_col = "feature_1",
    corr_feature2_col = "feature_2",
    corr_value_col = "rho"
) {
  
  assoc_sub <- assoc_long[
    assoc_long$family == family &
      assoc_long$pattern == pattern,
    ,
    drop = FALSE
  ]
  
  stopifnot(nrow(assoc_sub) > 0L)
  
  assoc_mat <- build_feature_pathway_matrix(
    assoc_long = assoc_sub,
    feature_col = feature_col,
    pathway_col = pathway_col,
    value_col = assoc_col
  )
  
  #======eliminar pathways sin variabilidad===========
  
  variable_pathways <- apply(
    assoc_mat,
    2,
    function(x) length(unique(x)) > 1L
  )
  
  assoc_mat <- assoc_mat[
    ,
    variable_pathways,
    drop = FALSE
  ]
  
  stopifnot(
    nrow(assoc_mat) >= 2L,
    ncol(assoc_mat) >= 2L
  )
  
  
  kept_features <- rownames(assoc_mat)
  kept_pathways <- colnames(assoc_mat)
  
  corr_sub <- corr_long[
    corr_long$family == family &
      corr_long$pattern == pattern &
      corr_long[[corr_feature1_col]] %in% kept_features &
      corr_long[[corr_feature2_col]] %in% kept_features,
    ,
    drop = FALSE
  ]
  
  feat_corr_mat <- build_symmetric_corr_matrix(
    corr_long = corr_sub,
    features = kept_features,
    feature1_col = corr_feature1_col,
    feature2_col = corr_feature2_col,
    value_col = corr_value_col
  )
  #======distancias desde correlaciones observadas===========
  
  missing_pairs <- sum(
    !is.finite(feat_corr_mat[upper.tri(feat_corr_mat)])
  )
  
  if (missing_pairs > 0L) {
    message(
      "AVISO: ", missing_pairs,
      " pares sin correlacion observada. ",
      "Se usara distancia neutra (=1) solo para el clustering piloto."
    )
  }
  
  feat_dist_matrix <- 1 - feat_corr_mat
  feat_dist_matrix[!is.finite(feat_dist_matrix)] <- 1
  diag(feat_dist_matrix) <- 0
  
  feat_dist <- stats::as.dist(feat_dist_matrix)
  
  feat_hc <- stats::hclust(
    feat_dist,
    method = PILOT_CLUSTER_METHOD
  )
  
  k_feat <- max(
    1L,
    min(
      as.integer(n_feature_modules),
      nrow(assoc_mat)
    )
  )
  
  feat_module <- cutree(feat_hc, k = k_feat)
  feat_module_df <- data.frame(
    feature = names(feat_module),
    feature_module = paste0("FM", unname(feat_module)),
    stringsAsFactors = FALSE
  )
  
  path_profile_corr <- suppressWarnings(
    cor(
      assoc_mat,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )
  
  if (is.null(dim(path_profile_corr))) {
    path_profile_corr <- matrix(
      1,
      nrow = 1,
      ncol = 1,
      dimnames = list(kept_pathways, kept_pathways)
    )
  }
  
  path_dist <- as.dist(1 - path_profile_corr)
  path_hc <- hclust(path_dist, method = PILOT_CLUSTER_METHOD)
  
  k_path <- max(
    1L,
    min(
      as.integer(n_pathway_modules),
      ncol(assoc_mat)
    )
  )
  
  path_module <- cutree(path_hc, k = k_path)
  path_module_df <- data.frame(
    pathway = names(path_module),
    pathway_module = paste0("PM", unname(path_module)),
    stringsAsFactors = FALSE
  )
  
  assoc_cells <- expand.grid(
    feature = rownames(assoc_mat),
    pathway = colnames(assoc_mat),
    stringsAsFactors = FALSE
  )
  
  assoc_cells$association_percentile <- as.vector(assoc_mat)
  #======presencia real de enlace Q5===========
  
  assoc_cells$has_link <- paste(
    assoc_cells$feature,
    assoc_cells$pathway,
    sep = "\r"
  ) %in% paste(
    as.character(assoc_sub[[feature_col]]),
    as.character(assoc_sub[[pathway_col]]),
    sep = "\r"
  )
  assoc_cells <- merge(
    assoc_cells,
    feat_module_df,
    by = "feature",
    all.x = TRUE
  )
  
  assoc_cells <- merge(
    assoc_cells,
    path_module_df,
    by = "pathway",
    all.x = TRUE
  )
  
  module_heat <- aggregate(
    association_percentile ~ feature_module + pathway_module,
    data = assoc_cells,
    FUN = mean
  )
  
  link_count <- stats::aggregate(
    has_link ~ feature_module + pathway_module,
    data = assoc_cells,
    FUN = sum
  )
  
  names(link_count)[3L] <- "n_links"
  
  feature_count <- aggregate(
    feature ~ feature_module,
    data = feat_module_df,
    FUN = function(x) length(unique(x))
  )
  names(feature_count)[2] <- "n_features"
  
  pathway_count <- aggregate(
    pathway ~ pathway_module,
    data = path_module_df,
    FUN = function(x) length(unique(x))
  )
  names(pathway_count)[2] <- "n_pathways"
  
  module_heat <- merge(
    module_heat,
    link_count,
    by = c("feature_module", "pathway_module"),
    all.x = TRUE
  )
  
  module_heat <- merge(
    module_heat,
    feature_count,
    by = "feature_module",
    all.x = TRUE
  )
  
  module_heat <- merge(
    module_heat,
    pathway_count,
    by = "pathway_module",
    all.x = TRUE
  )
  
  feat_levels <- unique(feat_module_df$feature_module[order(feat_module_df$feature_module)])
  path_levels <- unique(path_module_df$pathway_module[order(path_module_df$pathway_module)])
  
  module_heat$feature_module_label <- paste0(
    module_heat$feature_module,
    " (nF=",
    module_heat$n_features,
    ")"
  )
  
  module_heat$pathway_module_label <- paste0(
    module_heat$pathway_module,
    " (nP=",
    module_heat$n_pathways,
    ")"
  )
  
  feat_label_levels <- unique(
    paste0(feature_count$feature_module, " (nF=", feature_count$n_features, ")")
  )
  feat_label_levels <- feat_label_levels[order(feature_count$feature_module)]
  
  path_label_levels <- unique(
    paste0(pathway_count$pathway_module, " (nP=", pathway_count$n_pathways, ")")
  )
  path_label_levels <- path_label_levels[order(pathway_count$pathway_module)]
  
  module_heat$feature_module_label <- factor(
    module_heat$feature_module_label,
    levels = rev(feat_label_levels)
  )
  
  module_heat$pathway_module_label <- factor(
    module_heat$pathway_module_label,
    levels = path_label_levels
  )
  
  list(
    assoc_mat = assoc_mat,
    feat_corr_mat = feat_corr_mat,
    path_profile_corr = path_profile_corr,
    feat_module_df = feat_module_df,
    path_module_df = path_module_df,
    module_heat = module_heat,
    assoc_cells = assoc_cells
  )
}

#======funcion plot modulo x modulo===========
plot_feature_function_module_heatmap <- function(
    module_obj,
    family,
    pattern
) {
  
  ggplot2::ggplot(
    module_obj$module_heat,
    ggplot2::aes(
      x = pathway_module_label,
      y = feature_module_label,
      fill = association_percentile
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(
          "n=", n_links,
          "\n",
          sprintf("%.2f", association_percentile)
        )
      ),
      size = HM_CELL_TEXT_SIZE,
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_gradient(
      low = "grey95",
      high = "black",
      limits = c(0, 1),
      name = "Mean\nassociation\npercentile"
    ) +
    ggplot2::labs(
      title = paste0(
        family,
        " — ",
        pattern,
        ": feature modules × functional modules"
      ),
      subtitle = paste0(
        "Feature modules = clustering by OBSERVED feature-feature correlation. ",
        "Pathway modules = clustering by pathway-association profiles.\n",
        "Cell text: n = non-zero feature-pathway links; second line = mean association percentile."
      ),
      x = "Functional module",
      y = "Feature module"
    ) +
    ggplot2::theme_minimal(base_size = HM_BASE_SIZE) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = HM_X_ANGLE, hjust = 1, vjust = 1,
                                          size = HM_AXIS_X_SIZE),
      axis.text.y = ggplot2::element_text(size = HM_AXIS_Y_SIZE),
      plot.title = ggplot2::element_text(size = HM_TITLE_SIZE, face = "bold"),
      legend.title = ggplot2::element_text(size = HM_LEGEND_TITLE_SIZE),
      legend.text = ggplot2::element_text(size = HM_LEGEND_TEXT_SIZE)
    )
}

#======EJECUCION===========

SELECTION_RUN <- canonical(
  SELECTION_RUN
)

prepare_plot_features <- function(
    feature_table,
    sharing_table,
    subset_type = c(
      "all",
      "specific",
      "shared"
    )
) {
  #======preparar datos para barplots===========
  
  subset_type <- match.arg(
    subset_type
  )
  
  x <- feature_table
  
  idx <- match(
    x$feature_id,
    sharing_table$feature_id
  )
  
  assert(
    all(!is.na(idx)),
    "Hay features sin correspondencia en FEATURES_FINAL$SHARING."
  )
  
  x$n_questions <- sharing_table$n_questions[
    idx
  ]
  
  x$questions <- sharing_table$questions[
    idx
  ]
  
  #======conservar patrones FD NP SO===========
  
  x <- x[
    !is.na(
      x$framework_direction
    ) &
      x$framework_direction %in%
      DIRECTION_LEVELS,
    ,
    drop = FALSE
  ]
  
  #======filtrar especificas o compartidas===========
  
  if (subset_type == "specific") {
    x <- x[
      x$n_questions == 1L,
      ,
      drop = FALSE
    ]
  }
  
  if (subset_type == "shared") {
    x <- x[
      x$n_questions >= 2L,
      ,
      drop = FALSE
    ]
  }
  
  x$framework_direction <- factor(
    x$framework_direction,
    levels = DIRECTION_LEVELS
  )
  
  x$view_plot <- factor(
    x$view_code,
    levels = names(
      VIEW_LABELS
    ),
    labels = unname(
      VIEW_LABELS
    )
  )
  
  #======ordenar por evidence rank score===========
  x$framework_direction <- droplevels(
    x$framework_direction
  )
  
  x$view_plot <- droplevels(
    x$view_plot
  )
  x$feature_plot <- stats::reorder(
    x$feature_label,
    x$evidence_rank_score
  )
  
  x
}

make_feature_trellis <- function(
    feature_table,
    sharing_table,
    family,
    subset_type = c(
      "all",
      "specific",
      "shared"
    )
) {
  #======crear trellis horizontal FD NP SO===========
  
  assert(
    requireNamespace(
      "ggplot2",
      quietly = TRUE
    ),
    "Se necesita ggplot2."
  )
  
  subset_type <- match.arg(
    subset_type
  )
  
  x <- prepare_plot_features(
    feature_table = feature_table,
    sharing_table = sharing_table,
    subset_type = subset_type
  )
  
  subtitle_text <- switch(
    subset_type,
    all = "Todas las features finales",
    specific = "Features especificas de este escenario",
    shared = "Features compartidas con al menos otro escenario"
  )
  
  #======control de subconjunto vacio===========
  
  if (!nrow(x)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0,
          y = 0,
          label = paste(
            family,
            subtitle_text,
            "Sin features",
            sep = "\n"
          )
        ) +
        ggplot2::theme_void()
    )
  }
  
  #======barplot trellis===========
  
  ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = evidence_rank_score,
      y = feature_plot,
      fill = view_plot
    )
  ) +
    ggplot2::geom_col(
      width = 0.78
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(
        framework_direction
      ),
      scales = "free_y",
      space = "free_y",
      drop = TRUE
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        0,
        1
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +
    ggplot2::labs(
      title = family,
      subtitle = subtitle_text,
      x = "Evidence rank score",
      y = NULL,
      fill = "Omica"
    ) +
    ggplot2::theme_bw(base_size = BAR_BASE_SIZE) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = BAR_TEXT_SIZE),
      plot.title = ggplot2::element_text(size = BAR_TITLE_SIZE, face = "bold"),
      strip.text.y = ggplot2::element_text(size = BAR_STRIP_TEXT_SIZE, face = "bold"),
      legend.title = ggplot2::element_text(size = BAR_LEGEND_TITLE_SIZE),
      legend.text = ggplot2::element_text(size = BAR_LEGEND_TEXT_SIZE),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}



make_feature_group_barplot <- function(
    feature_table,
    sharing_table,
    family,
    direction,
    subset_type = c(
      "all",
      "shared"
    )
) {
  #======barplot horizontal de un solo grupo===========
  
  assert(
    requireNamespace(
      "ggplot2",
      quietly = TRUE
    ),
    "Se necesita ggplot2."
  )
  
  subset_type <- match.arg(
    subset_type
  )
  
  assert(
    direction %in% DIRECTION_LEVELS,
    paste(
      "Direccion no valida:",
      direction
    )
  )
  
  x <- prepare_plot_features(
    feature_table = feature_table,
    sharing_table = sharing_table,
    subset_type = subset_type
  )
  
  x <- x[
    as.character(
      x$framework_direction
    ) == direction,
    ,
    drop = FALSE
  ]
  
  subtitle_text <- switch(
    subset_type,
    all = "Todas las features finales",
    shared = "Features compartidas con al menos otro escenario"
  )
  
  if (!nrow(x)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0,
          y = 0,
          label = paste(
            family,
            direction,
            subtitle_text,
            "Sin features",
            sep = "\n"
          )
        ) +
        ggplot2::theme_void()
    )
  }
  
  #======orden dentro del grupo===========
  
  x <- x[
    order(
      x$evidence_rank_score,
      x$feature_label
    ),
    ,
    drop = FALSE
  ]
  
  x$feature_plot <- factor(
    x$feature_label,
    levels = unique(
      x$feature_label
    )
  )
  
  #======plot===========
  
  ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = evidence_rank_score,
      y = feature_plot,
      fill = view_plot
    )
  ) +
    ggplot2::geom_col(
      width = 0.78
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        0,
        1
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +
    ggplot2::labs(
      title = paste(
        family,
        direction,
        sep = " | "
      ),
      subtitle = subtitle_text,
      x = "Evidence rank score",
      y = NULL,
      fill = "Omica"
    ) +
    ggplot2::theme_bw(base_size = BAR_BASE_SIZE) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = BAR_TEXT_SIZE),
      plot.title = ggplot2::element_text(size = BAR_TITLE_SIZE, face = "bold"),
      legend.title = ggplot2::element_text(size = BAR_LEGEND_TITLE_SIZE),
      legend.text = ggplot2::element_text(size = BAR_LEGEND_TEXT_SIZE),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}


BIO_RUN <- find_latest_bio_run(
  SELECTION_RUN
)

PANEL10_RUN <- find_latest_model_run(
  bio_run = BIO_RUN,
  subdir = "FINAL_PREDICTIVE_PANEL",
  required_rel = file.path(
    "tables",
    "07_FINAL_PANEL_ALL_FAMILIES.csv"
  )
)

SGL16_RUN <- find_latest_model_run(
  bio_run = BIO_RUN,
  subdir = "FINAL_SPARSE_GROUP_LASSO_PANEL",
  required_rel = file.path(
    "tables",
    "08_FINAL_SPARSE_GROUP_PANEL_ALL_FAMILIES.csv"
  )
)

PLS17_RUN <- find_latest_model_run(
  bio_run = BIO_RUN,
  subdir = "FINAL_EVIDENCE_A_PLSDA_PCA",
  required_rel = file.path(
    "tables",
    "02_EVIDENCE_A_FEATURES_FOR_PLSDA_ALL_FAMILIES.csv"
  )
)

#======leer pesos nativos de los modelos MOFA finales===========

MOFA_IMPORTANCE_LIST <- lapply(
  FAMILIES,
  function(fam) {
    
    audit_path <- find_crossomic_input_audit(
      selection_run = SELECTION_RUN,
      family = fam
    )
    
    audit <- read_csv_checked(
      audit_path,
      required = c(
        "family",
        "model_path"
      )
    )
    
    assert(
      nrow(audit) >= 1L,
      paste(
        "00_INPUT_AUDIT.csv esta vacio para",
        fam
      )
    )
    
    hit <- audit[
      audit$family == fam,
      ,
      drop = FALSE
    ]
    
    assert(
      nrow(hit) >= 1L,
      paste(
        "No hay fila para",
        fam,
        "en",
        audit_path
      )
    )
    
    model_path <- as.character(
      hit$model_path[1L]
    )
    
    assert(
      !is.na(model_path) &&
        nzchar(model_path),
      paste(
        "model_path vacio para",
        fam
      )
    )
    
    extract_mofa_importance(
      model_path = model_path,
      family = fam
    )
  }
)

MOFA_IMPORTANCE <- do.call(
  rbind,
  MOFA_IMPORTANCE_LIST
)

rownames(
  MOFA_IMPORTANCE
) <- NULL

evidence <- read_csv_checked(
  file.path(
    PLS17_RUN,
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
    "framework_obs_evidence_grade"
  )
)

vip <- read_csv_checked(
  file.path(
    PLS17_RUN,
    "tables",
    "10_PLSDA_VIP_PRIORITY_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "predictor_key",
    "VIP_max_selected_components"
  )
)

pls_loadings <- read_csv_checked(
  file.path(
    PLS17_RUN,
    "tables",
    "11_PLSDA_LOADINGS_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "predictor_key",
    "component",
    "loading"
  )
)

pls_status <- read_csv_checked(
  file.path(
    PLS17_RUN,
    "tables",
    "14_MODEL_STATUS_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "chosen_ncomp"
  )
)

enet_panel <- read_csv_checked(
  file.path(
    PANEL10_RUN,
    "tables",
    "07_FINAL_PANEL_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "predictor_key"
  )
)

enet_stability <- read_csv_checked(
  file.path(
    PANEL10_RUN,
    "tables",
    "06_SELECTION_STABILITY_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "predictor_key",
    "stability_frequency"
  )
)

sgl_panel <- read_csv_checked(
  file.path(
    SGL16_RUN,
    "tables",
    "08_FINAL_SPARSE_GROUP_PANEL_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "predictor_key"
  )
)

sgl_stability <- read_csv_checked(
  file.path(
    SGL16_RUN,
    "tables",
    "06_FEATURE_STABILITY_AUDIT_ALL_FAMILIES.csv"
  ),
  required = c(
    "family",
    "predictor_key",
    "feature_stability"
  )
)

evidence$family <- as.character(
  evidence$family
)

evidence$predictor_key <- as.character(
  evidence$predictor_key
)

evidence$view_code <- tolower(
  as.character(
    evidence$view_code
  )
)

evidence$feature_original <- as.character(
  evidence$feature_original
)

evidence$feature_label <- as.character(
  evidence$feature_label
)

vip$VIP_max_selected_components <- suppressWarnings(
  as.numeric(
    vip$VIP_max_selected_components
  )
)

pls_loadings$component <- suppressWarnings(
  as.integer(
    pls_loadings$component
  )
)

pls_loadings$loading <- suppressWarnings(
  as.numeric(
    pls_loadings$loading
  )
)

enet_stability$stability_frequency <- suppressWarnings(
  as.numeric(
    enet_stability$stability_frequency
  )
)

sgl_stability$feature_stability <- suppressWarnings(
  as.numeric(
    sgl_stability$feature_stability
  )
)

PLS_DIRECTION <- build_pls_direction(
  pls17_run = PLS17_RUN,
  pls_loadings = pls_loadings,
  pls_status = pls_status
)

MASTER <- evidence

assert(
  all(
    MASTER$framework_rec_evidence_grade ==
      "A_strong"
  ),
  "Hay features que no son Evidence-A strong."
)

assert(
  all(
    MASTER$framework_obs_evidence_grade ==
      "A_strong"
  ),
  "Hay features observadas que no son Evidence-A strong."
)

assert(
  all(
    MASTER$framework_rec_bidirectional_pattern ==
      MASTER$framework_obs_bidirectional_pattern
  ),
  "Hay desacuerdo entre patron reconstruido y observado."
)

MASTER$framework_pattern <-
  MASTER$framework_rec_bidirectional_pattern

MASTER$framework_direction <-
  pattern_to_class(
    MASTER$framework_pattern
  )

MASTER <- merge(
  MASTER,
  MOFA_IMPORTANCE,
  by = c(
    "family",
    "view_code",
    "feature_original"
  ),
  all.x = TRUE,
  sort = FALSE
)

assert(
  all(
    is.finite(
      MASTER$MOFA_loading_norm
    )
  ),
  "Hay features Evidence-A sin peso MOFA. Revisar correspondencia feature/view."
)

MASTER <- merge(
  MASTER,
  vip[
    ,
    c(
      "family",
      "predictor_key",
      "VIP_max_selected_components"
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

MASTER <- merge(
  MASTER,
  PLS_DIRECTION,
  by = c(
    "family",
    "predictor_key"
  ),
  all.x = TRUE,
  sort = FALSE
)

#======voto PLS por importancia; la direccion la dicta solo el framework===========

MASTER$vote_PLSDA <-
  is.finite(
    MASTER$VIP_max_selected_components
  ) &
  MASTER$VIP_max_selected_components >=
  VIP_THRESHOLD

MASTER$PLS_framework_margin <- framework_pls_margin(
  direction = MASTER$framework_direction,
  align_np = MASTER$PLS_alignment_NP,
  align_fd = MASTER$PLS_alignment_FD,
  align_so = MASTER$PLS_alignment_SO
)

MASTER$vote_ElasticNet <-
  paste(
    MASTER$family,
    MASTER$predictor_key,
    sep = "\r"
  ) %in%
  paste(
    enet_panel$family,
    enet_panel$predictor_key,
    sep = "\r"
  )

MASTER$vote_SGL <-
  paste(
    MASTER$family,
    MASTER$predictor_key,
    sep = "\r"
  ) %in%
  paste(
    sgl_panel$family,
    sgl_panel$predictor_key,
    sep = "\r"
  )

enet_stability_small <- enet_stability[
  ,
  c(
    "family",
    "predictor_key",
    "stability_frequency"
  ),
  drop = FALSE
]

names(enet_stability_small)[3L] <-
  "ElasticNet_stability"

sgl_stability_small <- sgl_stability[
  ,
  c(
    "family",
    "predictor_key",
    "feature_stability"
  ),
  drop = FALSE
]

names(sgl_stability_small)[3L] <-
  "SGL_feature_stability"

MASTER <- merge(
  MASTER,
  enet_stability_small,
  by = c(
    "family",
    "predictor_key"
  ),
  all.x = TRUE,
  sort = FALSE
)

MASTER <- merge(
  MASTER,
  sgl_stability_small,
  by = c(
    "family",
    "predictor_key"
  ),
  all.x = TRUE,
  sort = FALSE
)

MASTER$n_model_votes <-
  as.integer(
    MASTER$vote_ElasticNet
  ) +
  as.integer(
    MASTER$vote_SGL
  ) +
  as.integer(
    MASTER$vote_PLSDA
  )

#======ranking multimetrica sin pesos arbitrarios===========

rank_group <- interaction(
  MASTER$family,
  MASTER$view_code,
  drop = TRUE
)

MASTER$MOFA_rank <- ave(
  MASTER$MOFA_loading_norm,
  rank_group,
  FUN = percentile_rank
)

MASTER$VIP_rank <- ave(
  MASTER$VIP_max_selected_components,
  rank_group,
  FUN = percentile_rank
)

MASTER$PLS_framework_margin_rank <- ave(
  MASTER$PLS_framework_margin,
  rank_group,
  FUN = percentile_rank
)

MASTER$PLS_framework_margin_rank[
  !is.finite(
    MASTER$PLS_framework_margin_rank
  )
] <- 0.5

MASTER$ElasticNet_stability_rank <- ave(
  ifelse(
    is.finite(
      MASTER$ElasticNet_stability
    ),
    MASTER$ElasticNet_stability,
    0
  ),
  rank_group,
  FUN = percentile_rank
)

MASTER$SGL_stability_rank <- ave(
  ifelse(
    is.finite(
      MASTER$SGL_feature_stability
    ),
    MASTER$SGL_feature_stability,
    0
  ),
  rank_group,
  FUN = percentile_rank
)

MASTER$evidence_rank_score <- rowMeans(
  cbind(
    MASTER$MOFA_rank,
    MASTER$VIP_rank,
    MASTER$PLS_framework_margin_rank,
    MASTER$ElasticNet_stability_rank,
    MASTER$SGL_stability_rank
  )
)

MASTER$support_tier <- ifelse(
  MASTER$n_model_votes == 3L,
  "3of3",
  ifelse(
    MASTER$n_model_votes == 2L,
    "2of3",
    paste0(
      MASTER$n_model_votes,
      "of3"
    )
  )
)

MASTER$feature_id <- paste(
  MASTER$view_code,
  MASTER$feature_original,
  sep = "::"
)

FINAL <- MASTER[
  MASTER$n_model_votes >=
    MIN_MODEL_VOTES,
  ,
  drop = FALSE
]

FINAL <- FINAL[
  order(
    FINAL$family,
    -FINAL$n_model_votes,
    -FINAL$evidence_rank_score,
    FINAL$view_code,
    FINAL$framework_pattern,
    FINAL$feature_label
  ),
  ,
  drop = FALSE
]

cols_final <- c(
  "family",
  "feature_id",
  "predictor_key",
  "view_code",
  "feature_original",
  "feature_label",
  "framework_pattern",
  "framework_direction",
  "support_tier",
  "n_model_votes",
  "vote_ElasticNet",
  "vote_SGL",
  "vote_PLSDA",
  "MOFA_loading_norm",
  "VIP_max_selected_components",
  "PLS_framework_margin",
  "ElasticNet_stability",
  "SGL_feature_stability",
  "MOFA_rank",
  "VIP_rank",
  "PLS_framework_margin_rank",
  "ElasticNet_stability_rank",
  "SGL_stability_rank",
  "evidence_rank_score"
)

FINAL <- FINAL[
  ,
  cols_final,
  drop = FALSE
]

P1 <- FINAL[
  FINAL$family == "P1",
  ,
  drop = FALSE
]

P2 <- FINAL[
  FINAL$family == "P2",
  ,
  drop = FALSE
]

P3 <- FINAL[
  FINAL$family == "P3",
  ,
  drop = FALSE
]

SHARING <- build_sharing_table(
  FINAL
)

FEATURES_FINAL <- list(
  P1 = P1,
  P2 = P2,
  P3 = P3,
  ALL = FINAL,
  SHARING = SHARING
)

if (
  interactive() &&
  isTRUE(OPEN_VIEWS)
) {
  View(
    FEATURES_FINAL$P1,
    title = "P1 - features finales"
  )
  
  View(
    FEATURES_FINAL$P2,
    title = "P2 - features finales"
  )
  
  View(
    FEATURES_FINAL$P3,
    title = "P3 - features finales"
  )
  
  View(
    FEATURES_FINAL$ALL,
    title = "P1 P2 P3 - todo junto"
  )
  
  View(
    FEATURES_FINAL$SHARING,
    title = "Comparacion P1 P2 P3"
  )
}

FEATURES_FINAL


#======PLOTS TODAS LAS FEATURES POR GRUPO===========

PLOTS_ALL <- list()

for (fam in FAMILIES) {
  
  PLOTS_ALL[[fam]] <- list()
  
  for (direction in DIRECTION_LEVELS) {
    
    PLOTS_ALL[[fam]][[direction]] <-
      make_feature_group_barplot(
        feature_table = FEATURES_FINAL[[fam]],
        sharing_table = FEATURES_FINAL$SHARING,
        family = fam,
        direction = direction,
        subset_type = "all"
      )
    
    print(
      PLOTS_ALL[[fam]][[direction]]
    )
  }
}

#======PLOTS FEATURES ESPECIFICAS===========

PLOTS_SPECIFIC <- list()

for (fam in FAMILIES) {
  
  PLOTS_SPECIFIC[[fam]] <- make_feature_trellis(
    feature_table = FEATURES_FINAL[[fam]],
    sharing_table = FEATURES_FINAL$SHARING,
    family = fam,
    subset_type = "specific"
  )
  
  print(
    PLOTS_SPECIFIC[[fam]]
  )
}



#======PLOTS FEATURES COMPARTIDAS POR GRUPO===========

PLOTS_SHARED <- list()

for (fam in FAMILIES) {
  
  PLOTS_SHARED[[fam]] <- list()
  
  for (direction in DIRECTION_LEVELS) {
    
    PLOTS_SHARED[[fam]][[direction]] <-
      make_feature_group_barplot(
        feature_table = FEATURES_FINAL[[fam]],
        sharing_table = FEATURES_FINAL$SHARING,
        family = fam,
        direction = direction,
        subset_type = "shared"
      )
    
    print(
      PLOTS_SHARED[[fam]][[direction]]
    )
  }
}


#======LISTA DE PLOTS===========

PLOTS_FINAL <- list(
  ALL = PLOTS_ALL,
  SPECIFIC = PLOTS_SPECIFIC,
  SHARED = PLOTS_SHARED
)

#======PARAMETROS PREVIEW RED TRIPARTITA===========

TRIP_PREVIEW_FAMILY <- "P1"
TRIP_PREVIEW_PATTERN <- "SO_like_high"
TRIP_PREVIEW_SCOPE <- "all"   # all / shared / specific

TRIP_KEEP_STRONGEST_PROGRAM_PER_FEATURE <- TRUE

TRIP_FEATURE_COLORS <- c(
  tx = "#B07AA1",
  pr = "#4E79A7",
  me = "#59A14F",
  cl = "#E15759"
)

TRIP_FEATURE_SIZE_RANGE <- c(3.2, 10.5)
TRIP_FP_EDGE_RANGE <- c(0.6, 2.4)
TRIP_PP_EDGE_RANGE <- c(1.0, 4.2)
TRIP_PROGRAM_NODE_SIZE <- 7.5
TRIP_PATHWAY_NODE_SIZE <- 4.5

TRIP_OPEN_PLOT <- TRUE


#======FUNCIONES RED TRIPARTITA===========

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

feature_key_trip <- function(view, original) {
  paste(
    tolower(as.character(view)),
    as.character(original),
    sep = "\r"
  )
}

program_key_trip <- function(family, side, program_id) {
  paste(
    family,
    side,
    program_id,
    sep = "\r"
  )
}

pathway_key_trip <- function(family, side, program_id, library, pathway) {
  paste(
    family,
    side,
    program_id,
    library,
    pathway,
    sep = "\r"
  )
}

rescale_to <- function(x, to = c(1, 6), default = mean(to)) {
  out <- rep(default, length(x))
  
  ok <- is.finite(x)
  
  if (!any(ok)) {
    return(out)
  }
  
  xr <- range(x[ok], na.rm = TRUE)
  
  if (xr[1] == xr[2]) {
    out[ok] <- default
    return(out)
  }
  
  out[ok] <- to[1] + (x[ok] - xr[1]) / (xr[2] - xr[1]) * (to[2] - to[1])
  
  out
}

normalize_join_key <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}
#======limpiar label visible de pathway===========

clean_pathway_label_trip <- function(
    x,
    width = 28
) {
  
  z <- as.character(
    x
  )
  
  # quitar "GO_BP |", "GO_CC |", "REACTOME |", etc.
  z <- sub(
    "^[^|]+\\|\\s*",
    "",
    z
  )
  
  # quitar prefijos repetidos dentro del nombre
  z <- sub(
    "^(GOBP|GO_BP|GOCC|GO_CC|GOMF|GO_MF|HALLMARK|REACTOME)_",
    "",
    z
  )
  
  # underscores -> espacios
  z <- gsub(
    "_+",
    " ",
    z
  )
  
  z <- trimws(
    z
  )
  
  # partir nombres largos en varias lineas
  vapply(
    z,
    function(s) {
      paste(
        strwrap(
          s,
          width = width
        ),
        collapse = "\n"
      )
    },
    character(1)
  )
}
extract_manifest_value_simple <- function(manifest_path, key) {
  assert(
    file.exists(manifest_path),
    paste("No existe RUN_MANIFEST:", manifest_path)
  )
  
  z <- utils::read.csv(
    manifest_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  nm <- tolower(names(z))
  
  key_col <- names(z)[match(TRUE, nm %in% c("key", "field", "item", "name"))]
  value_col <- names(z)[match(TRUE, nm %in% c("value", "path", "ruta"))]
  
  if (is.na(key_col) || is.na(value_col)) {
    if (ncol(z) >= 2L) {
      key_col <- names(z)[1L]
      value_col <- names(z)[2L]
    } else {
      stop("RUN_MANIFEST.csv no tiene formato interpretable.", call. = FALSE)
    }
  }
  
  hit <- z[z[[key_col]] == key, value_col, drop = TRUE]
  
  if (!length(hit)) {
    return(NA_character_)
  }
  
  as.character(hit[1L])
}

resolve_synthesis_run_trip <- function(panel10_run) {
  manifest_path <- file.path(
    panel10_run,
    "RUN_MANIFEST.csv"
  )
  
  p <- extract_manifest_value_simple(
    manifest_path,
    "BIO_SYNTHESIS_RUN"
  )
  
  assert(
    !is.na(p) && nzchar(p) && dir.exists(p),
    paste(
      "No se pudo resolver BIO_SYNTHESIS_RUN desde:",
      manifest_path
    )
  )
  
  canonical(p)
}

load_q5_catalog_trip <- function(panel10_run) {
  #======cargar catalogo congelado Q5===========
  synthesis_run <- resolve_synthesis_run_trip(
    panel10_run
  )
  
  q5_path <- file.path(
    synthesis_run,
    "plot_data",
    "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
    "00_INDEX",
    "FIG5_MASTER_PATHWAY_FEATURE_CATALOG.csv"
  )
  
  assert(
    file.exists(q5_path),
    paste("No existe FIG5_MASTER_PATHWAY_FEATURE_CATALOG.csv:", q5_path)
  )
  
  raw <- utils::read.csv(
    q5_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  req <- c(
    "family",
    "program_id",
    "library",
    "pathway",
    "phenotypic_side",
    "view_code",
    "feature_original"
  )
  
  miss <- setdiff(req, names(raw))
  
  assert(
    length(miss) == 0L,
    paste(
      "Faltan columnas en Q5 catalog:",
      paste(miss, collapse = ", ")
    )
  )
  
  out <- data.frame(
    family = as.character(raw$family),
    program_id = as.character(raw$program_id),
    program_short = as.character(raw$program_short %||% raw$program_id),
    library = as.character(raw$library),
    pathway = as.character(raw$pathway),
    pathway_label = paste(
      as.character(
        raw$library
      ),
      as.character(
        raw$pathway
      ),
      sep = " | "
    ),
    q5_phenotypic_side = as.character(raw$phenotypic_side),
    q5_view_code = tolower(as.character(raw$view_code)),
    q5_feature_original = as.character(raw$feature_original),
    q5_feature_label = as.character(raw$feature_label %||% raw$feature_original),
    tx_driver_weight = suppressWarnings(as.numeric(raw$tx_driver_weight %||% NA)),
    pathway_feature_cosine = suppressWarnings(as.numeric(raw$pathway_feature_cosine %||% NA)),
    stringsAsFactors = FALSE
  )
  
  out$program_label <- ifelse(
    !is.na(out$program_short) & nzchar(trimws(out$program_short)),
    out$program_short,
    out$program_id
  )
  
  out$feature_key <- feature_key_trip(
    out$q5_view_code,
    out$q5_feature_original
  )
  
  out$program_key <- program_key_trip(
    out$family,
    out$q5_phenotypic_side,
    out$program_id
  )
  
  out$pathway_key <- pathway_key_trip(
    out$family,
    out$q5_phenotypic_side,
    out$program_id,
    out$library,
    out$pathway
  )
  
  out$association_strength <- ifelse(
    out$q5_view_code == "tx" & is.finite(out$tx_driver_weight),
    abs(out$tx_driver_weight),
    abs(out$pathway_feature_cosine)
  )
  
  out$association_strength[!is.finite(out$association_strength)] <- 0
  
  out
}

load_optional_axis_maps_trip <- function() {
  #======cargar mapas auxiliares CL y PR si existen===========
  out <- list(
    cl = NULL,
    pr = NULL
  )
  
  cl_path <- file.path(
    getwd(),
    "reference_maps",
    "clinical_master_alias_map.tsv"
  )
  
  pr_path <- file.path(
    getwd(),
    "reference_maps",
    "protein_master_alias_map.tsv"
  )
  
  if (file.exists(cl_path)) {
    cl <- utils::read.delim(
      cl_path,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    key_col <- intersect(
      c("alias", "feature_label_es", "canonical_label_en"),
      names(cl)
    )[1L]
    
    if (!is.na(key_col) && "biological_axis" %in% names(cl)) {
      out$cl <- unique(
        data.frame(
          join_key = normalize_join_key(cl[[key_col]]),
          biological_axis = as.character(cl$biological_axis),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  if (file.exists(pr_path)) {
    pr <- utils::read.delim(
      pr_path,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    key_col <- intersect(
      c("alias", "protein_symbol", "gene_symbol", "feature_label"),
      names(pr)
    )[1L]
    
    if (!is.na(key_col) && "biological_axis" %in% names(pr)) {
      out$pr <- unique(
        data.frame(
          join_key = normalize_join_key(pr[[key_col]]),
          biological_axis = as.character(pr$biological_axis),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  out
}

annotate_axis_trip <- function(feature_table, axis_maps) {
  #======anotar biological_axis en CL y PR===========
  z <- feature_table
  z$biological_axis <- NA_character_
  
  for (vw in c("cl", "pr")) {
    map_now <- axis_maps[[vw]]
    
    if (is.null(map_now) || !nrow(map_now)) {
      next
    }
    
    idx <- which(z$view_code == vw)
    
    if (!length(idx)) {
      next
    }
    
    key1 <- normalize_join_key(z$feature_label[idx])
    key2 <- normalize_join_key(z$feature_original[idx])
    
    hit1 <- match(key1, map_now$join_key)
    hit2 <- match(key2, map_now$join_key)
    
    val1 <- map_now$biological_axis[hit1]
    val2 <- map_now$biological_axis[hit2]
    
    z$biological_axis[idx] <- ifelse(
      !is.na(val1) & nzchar(val1),
      val1,
      val2
    )
  }
  
  z
}

select_trip_features <- function(
    features_final,
    sharing_table,
    family,
    pattern,
    scope = c("all", "shared", "specific")
) {
  #======filtrar features segun familia patron y scope===========
  scope <- match.arg(scope)
  
  z <- features_final[[family]]
  
  z <- z[
    z$framework_pattern == pattern,
    ,
    drop = FALSE
  ]
  
  sh <- sharing_table[
    ,
    c(
      "feature_id",
      "n_questions",
      "questions"
    ),
    drop = FALSE
  ]
  
  z <- merge(
    z,
    sh,
    by = "feature_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  if (scope == "shared") {
    z <- z[
      z$n_questions > 1L,
      ,
      drop = FALSE
    ]
  }
  
  if (scope == "specific") {
    z <- z[
      z$n_questions == 1L &
        z$questions == family,
      ,
      drop = FALSE
    ]
  }
  
  z <- z[
    order(
      -z$evidence_rank_score,
      z$view_code,
      z$feature_label
    ),
    ,
    drop = FALSE
  ]
  
  rownames(z) <- NULL
  z
}

compute_program_y_trip <- function(
    feature_nodes,
    pathway_nodes,
    fp_edges,
    pp_edges,
    program_nodes
) {
  #======posicionar programas en y para reducir cruces===========
  raw_y <- rep(NA_real_, nrow(program_nodes))
  
  for (i in seq_len(nrow(program_nodes))) {
    pk <- program_nodes$program_key[i]
    
    fy <- feature_nodes$y[
      match(
        fp_edges$feature_id[
          fp_edges$program_key == pk
        ],
        feature_nodes$feature_id
      )
    ]
    
    py <- pathway_nodes$y[
      match(
        pp_edges$pathway_key[
          pp_edges$program_key == pk
        ],
        pathway_nodes$pathway_key
      )
    ]
    
    yy <- c(fy, py)
    yy <- yy[is.finite(yy)]
    
    raw_y[i] <- if (length(yy)) {
      mean(yy)
    } else {
      0
    }
  }
  
  ord <- order(raw_y, decreasing = TRUE)
  
  span_top <- max(
    c(feature_nodes$y, pathway_nodes$y),
    na.rm = TRUE
  )
  
  span_bottom <- min(
    c(feature_nodes$y, pathway_nodes$y),
    na.rm = TRUE
  )
  
  if (nrow(program_nodes) == 1L) {
    y_final <- mean(c(span_top, span_bottom))
  } else {
    y_final <- seq(
      from = span_top,
      to = span_bottom,
      length.out = nrow(program_nodes)
    )
  }
  
  out <- program_nodes
  out$y <- NA_real_
  out$y[ord] <- y_final
  
  out
}
#======limpiar nombres de pathways===========

clean_pathway_label <- function(x) {
  
  x <- as.character(x)
  
  x <- gsub("^GO_BP\\s*\\|\\s*", "", x)
  x <- gsub("^GO_CC\\s*\\|\\s*", "", x)
  x <- gsub("^HALLMARK\\s*\\|\\s*", "", x)
  x <- gsub("^REACTOME\\s*\\|\\s*", "", x)
  
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^GOCC_", "", x)
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^REACTOME_", "", x)
  
  x <- gsub("_", " ", x)
  x <- trimws(x)
  
  x
}

#======extraer libreria de pathway===========

extract_pathway_library <- function(x) {
  
  x <- as.character(x)
  
  out <- ifelse(
    grepl("^GO_BP\\s*\\|", x), "GO_BP",
    ifelse(
      grepl("^GO_CC\\s*\\|", x), "GO_CC",
      ifelse(
        grepl("^HALLMARK\\s*\\|", x), "HALLMARK",
        ifelse(
          grepl("^REACTOME\\s*\\|", x), "REACTOME",
          NA_character_
        )
      )
    )
  )
  
  out
}

#======partir texto en varias lineas===========

wrap_label_lines <- function(x, width = 26) {
  
  vapply(
    x,
    function(z) paste(strwrap(z, width = width), collapse = "\n"),
    character(1)
  )
}
build_tripartite_preview_data <- function(
    features_final,
    family,
    pattern,
    scope,
    q5_catalog,
    axis_maps,
    keep_strongest_program = TRUE
) {
  #======construir nodos y aristas===========
  feat <- select_trip_features(
    features_final = features_final,
    sharing_table = features_final$SHARING,
    family = family,
    pattern = pattern,
    scope = scope
  )
  
  assert(
    nrow(feat) > 0L,
    paste(
      "No hay features para",
      family,
      pattern,
      scope
    )
  )
  
  feat <- annotate_axis_trip(
    feat,
    axis_maps
  )
  
  feat$feature_key <- feature_key_trip(
    feat$view_code,
    feat$feature_original
  )
  
  side_now <- paste0(
    pattern_to_class(
      pattern
    ),
    "_side"
  )
  
  q5 <- q5_catalog[
    q5_catalog$family == family &
      q5_catalog$q5_phenotypic_side == side_now,
    ,
    drop = FALSE
  ]
  
  links <- merge(
    feat,
    q5,
    by = c("family", "feature_key"),
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )
  
  assert(
    nrow(links) > 0L,
    paste(
      "No hubo enlaces feature-program-pathway para",
      family,
      pattern,
      scope
    )
  )
  
  if (isTRUE(keep_strongest_program)) {
    prog_strength <- stats::aggregate(
      association_strength ~
        feature_id +
        program_key +
        program_id +
        program_label,
      data = links,
      FUN = max,
      na.rm = TRUE
    )
    
    keep_rows <- unlist(
      lapply(
        split(prog_strength, prog_strength$feature_id),
        function(d) {
          d <- d[
            order(
              -d$association_strength,
              d$program_id
            ),
            ,
            drop = FALSE
          ]
          rownames(d)[1L]
        }
      )
    )
    
    keep_prog <- prog_strength[
      rownames(prog_strength) %in% keep_rows,
      c("feature_id", "program_key"),
      drop = FALSE
    ]
    
    links <- merge(
      links,
      keep_prog,
      by = c("feature_id", "program_key"),
      all = FALSE,
      sort = FALSE
    )
  }
  
  fp_edges <- unique(
    links[
      ,
      c(
        "feature_id",
        "feature_label",
        "feature_original",
        "view_code",
        "support_tier",
        "n_model_votes",
        "evidence_rank_score",
        "biological_axis",
        "program_key",
        "program_id",
        "program_label",
        "association_strength"
      ),
      drop = FALSE
    ]
  )
  
  pp_edges <- stats::aggregate(
    feature_id ~
      program_key +
      program_id +
      program_label +
      pathway_key +
      pathway_label +
      library,
    data = links,
    FUN = function(x) length(unique(x))
  )
  
  names(pp_edges)[names(pp_edges) == "feature_id"] <- "n_retained_features"
  
  feature_nodes <- unique(
    fp_edges[
      ,
      c(
        "feature_id",
        "feature_label",
        "feature_original",
        "view_code",
        "support_tier",
        "n_model_votes",
        "evidence_rank_score",
        "biological_axis"
      ),
      drop = FALSE
    ]
  )
  
  feature_nodes <- feature_nodes[
    order(
      -feature_nodes$evidence_rank_score,
      feature_nodes$view_code,
      feature_nodes$feature_label
    ),
    ,
    drop = FALSE
  ]
  
  feature_nodes$y <- rev(
    seq_len(nrow(feature_nodes))
  )
  
  feature_nodes$x <- 0
  
  feature_nodes$node_size <- rescale_to(
    feature_nodes$evidence_rank_score,
    to = TRIP_FEATURE_SIZE_RANGE,
    default = mean(TRIP_FEATURE_SIZE_RANGE)
  )
  
  #======label limpio de feature===========
  
  feature_nodes$label <-
    feature_nodes$feature_label
  program_nodes <- unique(
    fp_edges[
      ,
      c(
        "program_key",
        "program_id",
        "program_label"
      ),
      drop = FALSE
    ]
  )
  program_nodes$display_label <- sub(
    "^.*__(P[0-9]+)$",
    "\\1",
    program_nodes$program_id
  )
  #======preparar nombres cortos de programas===========
  
  fp_edges$program_display <- sub(
    "^.*__(P[0-9]+)$",
    "\\1",
    fp_edges$program_id
  )
  
  pp_edges$program_display <- sub(
    "^.*__(P[0-9]+)$",
    "\\1",
    pp_edges$program_id
  )
  
  #======orden vertical de programas segun sus features===========
  
  program_nodes$feature_y_mean <- vapply(
    program_nodes$program_key,
    function(pk) {
      
      ids <- unique(
        fp_edges$feature_id[
          fp_edges$program_key == pk
        ]
      )
      
      yy <- feature_nodes$y[
        match(
          ids,
          feature_nodes$feature_id
        )
      ]
      
      mean(
        yy,
        na.rm = TRUE
      )
    },
    numeric(1)
  )
  
  program_nodes <- program_nodes[
    order(
      program_nodes$feature_y_mean,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  PROGRAM_DISPLAY_ORDER <-
    program_nodes$display_label
  
  #======crear nodos pathway===========
  
  pathway_nodes <- unique(
    pp_edges[
      ,
      c(
        "pathway_key",
        "pathway_label",
        "library"
      ),
      drop = FALSE
    ]
  )
  
  #======label corto para mostrar===========
  
  pathway_nodes$pathway_display <-
    clean_pathway_label_trip(
      pathway_nodes$pathway_label,
      width = 28
    )
  
  #======programa dominante de cada pathway===========
  
  pathway_dom <- pp_edges[
    ,
    c(
      "pathway_key",
      "program_display",
      "n_retained_features"
    ),
    drop = FALSE
  ]
  
  pathway_dom$program_rank <- match(
    pathway_dom$program_display,
    PROGRAM_DISPLAY_ORDER
  )
  
  pathway_dom <- pathway_dom[
    order(
      pathway_dom$pathway_key,
      -pathway_dom$n_retained_features,
      pathway_dom$program_rank
    ),
    ,
    drop = FALSE
  ]
  
  pathway_dom <- pathway_dom[
    !duplicated(
      pathway_dom$pathway_key
    ),
    ,
    drop = FALSE
  ]
  
  names(pathway_dom)[
    names(pathway_dom) ==
      "program_display"
  ] <- "dominant_program"
  
  names(pathway_dom)[
    names(pathway_dom) ==
      "n_retained_features"
  ] <- "dominant_support"
  
  pathway_nodes <- merge(
    pathway_nodes,
    pathway_dom[
      ,
      c(
        "pathway_key",
        "dominant_program",
        "dominant_support"
      ),
      drop = FALSE
    ],
    by = "pathway_key",
    all.x = TRUE,
    sort = FALSE
  )
  
  #======orden secundario por tipo de anotacion===========
  
  PATHWAY_LIBRARY_ORDER <- c(
    "HALLMARK",
    "REACTOME",
    "GO_BP",
    "GO_MF",
    "GO_CC"
  )
  
  pathway_nodes$program_order <- match(
    pathway_nodes$dominant_program,
    PROGRAM_DISPLAY_ORDER
  )
  
  pathway_nodes$library_order <- match(
    toupper(
      pathway_nodes$library
    ),
    PATHWAY_LIBRARY_ORDER
  )
  
  pathway_nodes$library_order[
    is.na(
      pathway_nodes$library_order
    )
  ] <- length(
    PATHWAY_LIBRARY_ORDER
  ) + 1L
  
  #======orden final: programa dominante -> libreria -> soporte -> nombre===========
  
  pathway_nodes <- pathway_nodes[
    order(
      pathway_nodes$program_order,
      pathway_nodes$library_order,
      -pathway_nodes$dominant_support,
      pathway_nodes$pathway_display
    ),
    ,
    drop = FALSE
  ]
  
  pathway_nodes$y <- rev(
    seq_len(
      nrow(
        pathway_nodes
      )
    )
  )
  
  pathway_nodes$x <- 2
  pathway_nodes$x <- 2
  
  #======posicionar programas respetando grupos de features===========
  
  program_nodes$y <- seq(
    from = max(
      c(
        feature_nodes$y,
        pathway_nodes$y
      ),
      na.rm = TRUE
    ),
    to = min(
      c(
        feature_nodes$y,
        pathway_nodes$y
      ),
      na.rm = TRUE
    ),
    length.out = nrow(
      program_nodes
    )
  )
  
  program_nodes$x <- 1
  program_nodes$x <- 1
  
  fp_edges <- merge(
    fp_edges,
    feature_nodes[
      ,
      c(
        "feature_id",
        "x",
        "y",
        "view_code"
      ),
      drop = FALSE
    ],
    by = c("feature_id", "view_code"),
    all.x = TRUE,
    sort = FALSE
  )
  
  names(fp_edges)[names(fp_edges) == "x"] <- "x_feature"
  names(fp_edges)[names(fp_edges) == "y"] <- "y_feature"
  
  fp_edges <- merge(
    fp_edges,
    program_nodes[
      ,
      c(
        "program_key",
        "x",
        "y"
      ),
      drop = FALSE
    ],
    by = "program_key",
    all.x = TRUE,
    sort = FALSE
  )
  
  names(fp_edges)[names(fp_edges) == "x"] <- "x_program"
  names(fp_edges)[names(fp_edges) == "y"] <- "y_program"
  
  fp_edges$edge_size <- rescale_to(
    fp_edges$association_strength,
    to = TRIP_FP_EDGE_RANGE,
    default = mean(TRIP_FP_EDGE_RANGE)
  )
  
  pp_edges <- merge(
    pp_edges,
    program_nodes[
      ,
      c(
        "program_key",
        "x",
        "y"
      ),
      drop = FALSE
    ],
    by = "program_key",
    all.x = TRUE,
    sort = FALSE
  )
  
  names(pp_edges)[names(pp_edges) == "x"] <- "x_program"
  names(pp_edges)[names(pp_edges) == "y"] <- "y_program"
  
  pp_edges <- merge(
    pp_edges,
    pathway_nodes[
      ,
      c(
        "pathway_key",
        "x",
        "y"
      ),
      drop = FALSE
    ],
    by = "pathway_key",
    all.x = TRUE,
    sort = FALSE
  )
  
  names(pp_edges)[names(pp_edges) == "x"] <- "x_pathway"
  names(pp_edges)[names(pp_edges) == "y"] <- "y_pathway"
  
  pp_edges$edge_size <- rescale_to(
    pp_edges$n_retained_features,
    to = TRIP_PP_EDGE_RANGE,
    default = mean(TRIP_PP_EDGE_RANGE)
  )
  
  list(
    feature_nodes = feature_nodes,
    program_nodes = program_nodes,
    pathway_nodes = pathway_nodes,
    fp_edges = fp_edges,
    pp_edges = pp_edges
  )
}

plot_tripartite_preview <- function(
    trip_data,
    family,
    pattern,
    scope
) {
  #======plot tripartito preview en una sola pagina===========
  assert(
    requireNamespace("ggplot2", quietly = TRUE),
    "Se necesita ggplot2."
  )
  
  fz <- trip_data$feature_nodes
  pz <- trip_data$program_nodes
  wz <- trip_data$pathway_nodes
  
  efp <- trip_data$fp_edges
  epp <- trip_data$pp_edges
  
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = efp,
      ggplot2::aes(
        x = x_feature,
        y = y_feature,
        xend = x_program,
        yend = y_program,
        linewidth = edge_size,
        colour = program_display
      ),
      linetype = "solid",
      alpha = 0.78,
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = epp,
      ggplot2::aes(
        x = x_program,
        y = y_program,
        xend = x_pathway,
        yend = y_pathway,
        linewidth = edge_size,
        colour = program_display
      ),
      alpha = 0.42,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = fz,
      ggplot2::aes(
        x = x,
        y = y,
        size = node_size,
        fill = view_code
      ),
      shape = 21,
      color = "black",
      stroke = 0.3
    ) +
    ggplot2::geom_point(
      data = pz,
      ggplot2::aes(
        x = x,
        y = y
      ),
      shape = 24,
      size = TRIP_PROGRAM_NODE_SIZE,
      fill = "grey75",
      color = "black",
      stroke = 0.4
    ) +
    ggplot2::geom_point(
      data = wz,
      ggplot2::aes(
        x = x,
        y = y,
        shape = library
      ),
      size = TRIP_PATHWAY_NODE_SIZE,
      colour = "black",
      stroke = 0.3
    ) +
    ggplot2::geom_text(
      data = fz,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 1,
      nudge_x = -0.03,
      size = TRIP_FEATURE_LABEL_SIZE,
      lineheight = 0.92
    ) +
    ggplot2::geom_text(
      data = pz,
      ggplot2::aes(
        x = x,
        y = y,
        label = display_label
      ),
      hjust = 0.5,
      nudge_y = 0.20,
      size = TRIP_PROGRAM_LABEL_SIZE,
      fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = wz,
      ggplot2::aes(
        x = x,
        y = y,
        label = pathway_display
      ),
      hjust = 0,
      nudge_x = 0.03,
      size = TRIP_PATHWAY_LABEL_SIZE,
      lineheight = 0.90
    ) +
    ggplot2::scale_fill_manual(
      values = TRIP_FEATURE_COLORS,
      labels = c(
        tx = "TX",
        pr = "PR",
        me = "ME",
        cl = "CL"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = EDGE_PROGRAM_COLORS,
      name = "Q5 program"
    ) +
    ggplot2::scale_shape_discrete(
      name = "Annotation source"
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_x_continuous(
      breaks = c(0, 1, 2),
      labels = c(
        "Retained feature",
        "Q5 program",
        "Exact pathway"
      ),
      expand = ggplot2::expansion(mult = c(0.18, 0.18))
    ) +
    ggplot2::coord_cartesian(
      clip = "off"
    ) +
    ggplot2::labs(
      title = paste0(
        family,
        " — ",
        pattern,
        ": tripartite biological story"
      ),
      subtitle = paste0(
        "Scope = ",
        scope,
        ". Feature -> strongest Q5 program -> associated biological annotation."
      ),
      x = NULL,
      y = NULL,
      fill = "Feature view"
    ) +
    ggplot2::theme_minimal(
      base_size = TRIP_BASE_SIZE
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = TRIP_LEGEND_TITLE_SIZE),
      legend.text = ggplot2::element_text(size = TRIP_LEGEND_TEXT_SIZE),
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      plot.margin = grid::unit(
        TRIP_PLOT_MARGIN,
        "lines"
      )
    )
}


#======EJECUCION PREVIEW RED TRIPARTITA===========

AXIS_MAPS_TRIP <- load_optional_axis_maps_trip()

Q5_CATALOG_TRIP <- load_q5_catalog_trip(
  panel10_run = PANEL10_RUN
)

TRIP_DATA_PREVIEW <- build_tripartite_preview_data(
  features_final = FEATURES_FINAL,
  family = TRIP_PREVIEW_FAMILY,
  pattern = TRIP_PREVIEW_PATTERN,
  scope = TRIP_PREVIEW_SCOPE,
  q5_catalog = Q5_CATALOG_TRIP,
  axis_maps = AXIS_MAPS_TRIP,
  keep_strongest_program = TRIP_KEEP_STRONGEST_PROGRAM_PER_FEATURE
)

PLOT_TRIPARTITE_PREVIEW <- plot_tripartite_preview(
  trip_data = TRIP_DATA_PREVIEW,
  family = TRIP_PREVIEW_FAMILY,
  pattern = TRIP_PREVIEW_PATTERN,
  scope = TRIP_PREVIEW_SCOPE
)

if (isTRUE(TRIP_OPEN_PLOT)) {
  print(PLOT_TRIPARTITE_PREVIEW)
}
#======TODAS LAS REDES TRIPARTITAS===========

TRIP_FAMILIES <- c(
  "P1",
  "P2",
  "P3"
)

TRIP_PATTERNS <- c(
  "FD_like_high",
  "NP_like_high",
  "SO_like_high"
)

TRIP_SCOPES <- c(
  "all",
  "shared",
  "specific"
)


#======contenedores===========

TRIPARTITE_DATA <- list()

TRIPARTITE_PLOTS <- list(
  ALL = list(),
  SHARED = list(),
  SPECIFIC = list()
)


#======generar todas las combinaciones===========

for (fam in TRIP_FAMILIES) {
  
  TRIPARTITE_DATA[[fam]] <- list()
  
  TRIPARTITE_PLOTS$ALL[[fam]] <- list()
  TRIPARTITE_PLOTS$SHARED[[fam]] <- list()
  TRIPARTITE_PLOTS$SPECIFIC[[fam]] <- list()
  
  for (pattern in TRIP_PATTERNS) {
    
    TRIPARTITE_DATA[[fam]][[pattern]] <- list()
    
    for (scope_now in TRIP_SCOPES) {
      
      cat(
        "\n============================================================\n"
      )
      
      cat(
        "TRIPARTITE:",
        fam,
        "|",
        pattern,
        "|",
        scope_now,
        "\n"
      )
      
      cat(
        "============================================================\n"
      )
      
      
      #======construir datos; si no existen features se salta===========
      
      trip_now <- tryCatch(
        {
          build_tripartite_preview_data(
            features_final = FEATURES_FINAL,
            family = fam,
            pattern = pattern,
            scope = scope_now,
            q5_catalog = Q5_CATALOG_TRIP,
            axis_maps = AXIS_MAPS_TRIP,
            keep_strongest_program =
              TRIP_KEEP_STRONGEST_PROGRAM_PER_FEATURE
          )
        },
        error = function(e) {
          
          message(
            "SKIP: ",
            fam,
            " | ",
            pattern,
            " | ",
            scope_now,
            " -> ",
            conditionMessage(e)
          )
          
          NULL
        }
      )
      
      
      #======si esta vacio continuar===========
      
      if (is.null(trip_now)) {
        next
      }
      
      
      #======guardar datos en memoria===========
      
      TRIPARTITE_DATA[[fam]][[pattern]][[scope_now]] <-
        trip_now
      
      
      #======crear plot===========
      
      plot_now <- plot_tripartite_preview(
        trip_data = trip_now,
        family = fam,
        pattern = pattern,
        scope = scope_now
      )
      
      
      #======guardar plot en memoria===========
      
      if (scope_now == "all") {
        
        TRIPARTITE_PLOTS$ALL[[fam]][[pattern]] <-
          plot_now
      }
      
      if (scope_now == "shared") {
        
        TRIPARTITE_PLOTS$SHARED[[fam]][[pattern]] <-
          plot_now
      }
      
      if (scope_now == "specific") {
        
        TRIPARTITE_PLOTS$SPECIFIC[[fam]][[pattern]] <-
          plot_now
      }
      
      
      #======mostrar en RStudio===========
      
      print(
        plot_now
      )
    }
  }
}

#======OBJETO FINAL DE REDES===========

NETWORKS_FINAL <- list(
  DATA = TRIPARTITE_DATA,
  PLOTS = TRIPARTITE_PLOTS
)

NETWORKS_FINAL


#======TODOS LOS HEATMAPS P1 P2 P3===========

HEATMAP_OPEN_TABLES <- FALSE
HEATMAP_RESULTS <- list()

for (fam_now in FAMILIES) {
  
  HEATMAP_RESULTS[[fam_now]] <- list()
  
  for (pattern_now in c(
    "FD_like_high",
    "NP_like_high",
    "SO_like_high"
  )) {
    
    message(
      "\nHEATMAPS: ", fam_now, " | ", pattern_now
    )
    
    HEATMAP_RESULTS[[fam_now]][pattern_now] <- list(
      tryCatch(
        local({
          
          HEAT_FAMILY <- fam_now
          HEAT_PATTERN <- pattern_now
          HEAT_SCOPE <- "all"
          
          
          #======funcion percentil robusto para links===========
          
          percentile_link <- function(x) {
            
            out <- rep(
              NA_real_,
              length(x)
            )
            
            ok <- is.finite(
              x
            )
            
            if (!any(ok)) {
              return(out)
            }
            
            z <- x[ok]
            
            if (
              length(z) == 1L ||
              length(unique(z)) == 1L
            ) {
              out[ok] <- 0.5
              return(out)
            }
            
            out[ok] <- (
              rank(
                z,
                ties.method = "average"
              ) - 1
            ) / (
              length(z) - 1
            )
            
            out
          }
          
          
          #======usar la misma seleccion feature -> programa de la red===========
          
          HEAT_TRIP_DATA <- build_tripartite_preview_data(
            features_final = FEATURES_FINAL,
            family = HEAT_FAMILY,
            pattern = HEAT_PATTERN,
            scope = HEAT_SCOPE,
            q5_catalog = Q5_CATALOG_TRIP,
            axis_maps = AXIS_MAPS_TRIP,
            keep_strongest_program = TRUE
          )
          
          
          #======programa elegido para cada feature===========
          
          HEAT_SELECTED_PROGRAM <- unique(
            HEAT_TRIP_DATA$fp_edges[
              ,
              c(
                "feature_id",
                "feature_label",
                "feature_original",
                "view_code",
                "evidence_rank_score",
                "program_key",
                "program_id"
              ),
              drop = FALSE
            ]
          )
          
          
          #======preparar Q5 para el mismo lado fenotipico===========
          
          HEAT_SIDE <- paste0(
            pattern_to_class(
              HEAT_PATTERN
            ),
            "_side"
          )
          
          HEAT_Q5 <- Q5_CATALOG_TRIP[
            Q5_CATALOG_TRIP$family ==
              HEAT_FAMILY &
              Q5_CATALOG_TRIP$q5_phenotypic_side ==
              HEAT_SIDE,
            ,
            drop = FALSE
          ]
          
          
          #======crear feature_id compatible===========
          
          HEAT_Q5$feature_id <- paste(
            HEAT_Q5$q5_view_code,
            HEAT_Q5$q5_feature_original,
            sep = "::"
          )
          
          
          #======quedarnos solo con links del programa elegido===========
          
          HEAT_LINKS <- merge(
            HEAT_SELECTED_PROGRAM,
            HEAT_Q5[
              ,
              c(
                "feature_id",
                "program_key",
                "library",
                "pathway",
                "pathway_label",
                "association_strength"
              ),
              drop = FALSE
            ],
            by = c(
              "feature_id",
              "program_key"
            ),
            all = FALSE,
            sort = FALSE
          )
          
          
          #======eliminar duplicados feature-pathway===========
          
          HEAT_LINKS <- stats::aggregate(
            association_strength ~
              feature_id +
              feature_label +
              feature_original +
              view_code +
              evidence_rank_score +
              program_id +
              library +
              pathway +
              pathway_label,
            data = HEAT_LINKS,
            FUN = max,
            na.rm = TRUE
          )
          
          
          #======normalizar fuerza de asociacion dentro de cada omica===========
          
          HEAT_LINKS$association_percentile <- ave(
            HEAT_LINKS$association_strength,
            HEAT_LINKS$view_code,
            FUN = percentile_link
          )
          
          
          #======nombre corto del programa===========
          
          HEAT_LINKS$program_display <- sub(
            "^.*__(P[0-9]+)$",
            "\\1",
            HEAT_LINKS$program_id
          )
          
          
          #======nombre limpio pathway===========
          
          HEAT_LINKS$pathway_display <-
            clean_pathway_label_trip(
              HEAT_LINKS$pathway_label,
              width = 22
            )
          
          
          #======mantener fuente de anotacion visible pero compacta===========
          
          HEAT_LINKS$pathway_display <- paste0(
            HEAT_LINKS$library,
            "\n",
            HEAT_LINKS$pathway_display
          )
          
          
          #======orden de features por ranking final===========
          
          FEATURE_ORDER <- unique(
            HEAT_SELECTED_PROGRAM[
              order(
                HEAT_SELECTED_PROGRAM$evidence_rank_score,
                HEAT_SELECTED_PROGRAM$feature_label
              ),
              "feature_label"
            ]
          )
          
          HEAT_LINKS$feature_label <- factor(
            HEAT_LINKS$feature_label,
            levels = FEATURE_ORDER
          )
          
          
          #======ordenar pathways por programa y libreria===========
          
          PATH_LIBRARY_ORDER <- c(
            "HALLMARK",
            "REACTOME",
            "GO_BP",
            "GO_MF",
            "GO_CC"
          )
          
          HEAT_LINKS$library_order <- match(
            HEAT_LINKS$library,
            PATH_LIBRARY_ORDER
          )
          
          HEAT_LINKS$library_order[
            is.na(
              HEAT_LINKS$library_order
            )
          ] <- length(
            PATH_LIBRARY_ORDER
          ) + 1L
          
          
          PATH_ORDER_TABLE <- unique(
            HEAT_LINKS[
              ,
              c(
                "program_display",
                "library_order",
                "pathway_display"
              ),
              drop = FALSE
            ]
          )
          
          PATH_ORDER_TABLE <- PATH_ORDER_TABLE[
            order(
              PATH_ORDER_TABLE$program_display,
              PATH_ORDER_TABLE$library_order,
              PATH_ORDER_TABLE$pathway_display
            ),
            ,
            drop = FALSE
          ]
          
          PATH_ORDER <- unique(
            PATH_ORDER_TABLE$pathway_display
          )
          
          HEAT_LINKS$pathway_display <- factor(
            HEAT_LINKS$pathway_display,
            levels = PATH_ORDER
          )
          
          
          #======crear matriz completa feature x pathway===========
          
          HEAT_GRID <- expand.grid(
            feature_label = levels(
              HEAT_LINKS$feature_label
            ),
            pathway_display = levels(
              HEAT_LINKS$pathway_display
            ),
            stringsAsFactors = FALSE
          )
          
          HEAT_GRID <- merge(
            HEAT_GRID,
            HEAT_LINKS[
              ,
              c(
                "feature_label",
                "pathway_display",
                "association_percentile",
                "view_code",
                "program_display"
              ),
              drop = FALSE
            ],
            by = c(
              "feature_label",
              "pathway_display"
            ),
            all.x = TRUE,
            sort = FALSE
          )
          
          
          #======restaurar factores===========
          
          HEAT_GRID$feature_label <- factor(
            HEAT_GRID$feature_label,
            levels = FEATURE_ORDER
          )
          
          HEAT_GRID$pathway_display <- factor(
            HEAT_GRID$pathway_display,
            levels = PATH_ORDER
          )
          
          
          #======PLOT HEATMAP===========
          
          PLOT_HEATMAP_PILOT <- ggplot2::ggplot(
            HEAT_GRID,
            ggplot2::aes(
              x = pathway_display,
              y = feature_label,
              fill = association_percentile
            )
          ) +
            
            ggplot2::geom_tile(
              colour = "white",
              linewidth = 0.25
            ) +
            
            ggplot2::scale_fill_gradient(
              low = "white",
              high = "black",
              limits = c(
                0,
                1
              ),
              na.value = "grey95",
              name = "Q5 association\npercentile"
            ) +
            
            ggplot2::labs(
              title = paste0(
                HEAT_FAMILY,
                " — ",
                HEAT_PATTERN,
                ": feature–pathway association map"
              ),
              subtitle = paste0(
                "Rows = retained features.\n",
                "Columns = pathways associated through each feature's strongest Q5 program.\n",
                "Color = within-view percentile of Q5 association strength."
              ),
              x = NULL,
              y = NULL
            ) +
            
            ggplot2::theme_minimal(base_size = HM_BASE_SIZE) +
            ggplot2::theme(
              panel.grid = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_text(
                angle = HM_X_ANGLE, hjust = 1, vjust = 1, size = 21
              ),
              axis.text.y = ggplot2::element_text(size = 23),
              legend.position = "right",
              legend.title = ggplot2::element_text(size = HM_LEGEND_TITLE_SIZE),
              legend.text = ggplot2::element_text(size = HM_LEGEND_TEXT_SIZE),
              plot.title = ggplot2::element_text(face = "bold", size = HM_TITLE_SIZE),
              plot.subtitle = ggplot2::element_text(size = 15, lineheight = 1.1,
                                                    margin = ggplot2::margin(b = 12)),
              plot.margin = ggplot2::margin(20, 25, 25, 20)
            )
          
          
          #======mostrar piloto===========
          
          print(
            PLOT_HEATMAP_PILOT
          )
          
          #======PILOTO CORRELACION FEATURE-FEATURE OBSERVED vs LATENT===========
          
          CORR_FAMILY <- HEAT_FAMILY
          CORR_PATTERN <- HEAT_PATTERN
          CORR_SCOPE <- HEAT_SCOPE
          
          # GLOBAL = pooled descriptivo Q5
          CORR_SAMPLE_SCOPE <- "GLOBAL"
          
          
          #======usar la misma red tripartita ya construida===========
          
          if (
            exists("TRIPARTITE_DATA") &&
            !is.null(
              TRIPARTITE_DATA[[CORR_FAMILY]][[CORR_PATTERN]][[CORR_SCOPE]]
            )
          ) {
            
            CORR_TRIP <-
              TRIPARTITE_DATA[[CORR_FAMILY]][[CORR_PATTERN]][[CORR_SCOPE]]
            
          } else {
            
            CORR_TRIP <- build_tripartite_preview_data(
              features_final = FEATURES_FINAL,
              family = CORR_FAMILY,
              pattern = CORR_PATTERN,
              scope = CORR_SCOPE,
              q5_catalog = Q5_CATALOG_TRIP,
              axis_maps = AXIS_MAPS_TRIP,
              keep_strongest_program = TRUE
            )
          }
          
          
          #======features retenidas del plot===========
          
          CORR_FEATURES <- unique(
            CORR_TRIP$feature_nodes[
              ,
              c(
                "feature_id",
                "feature_label",
                "feature_original",
                "view_code",
                "evidence_rank_score"
              ),
              drop = FALSE
            ]
          )
          
          CORR_FEATURES$matrix_key <- paste(
            CORR_FEATURES$view_code,
            CORR_FEATURES$feature_original,
            sep = "\r"
          )
          
          CORR_FEATURES <- CORR_FEATURES[
            order(
              -CORR_FEATURES$evidence_rank_score,
              CORR_FEATURES$feature_label
            ),
            ,
            drop = FALSE
          ]
          
          
          
          #======localizar matrices FIG5B ya calculadas===========
          
          CORR_SYNTHESIS_RUN <-
            resolve_synthesis_run_trip(
              PANEL10_RUN
            )
          
          CORR_ROOT <- file.path(
            CORR_SYNTHESIS_RUN,
            "plot_data",
            "FIG5_MULTIMODAL_PATHWAY_ARCHITECTURE",
            "01_ALL",
            CORR_FAMILY
          )
          
          
          OBS_FILES <- list.files(
            CORR_ROOT,
            pattern = paste0(
              "^FIG5B_OBSERVED_",
              CORR_SAMPLE_SCOPE,
              "_FEATURE_FEATURE_SPEARMAN_RHO_MATRIX\\.csv$"
            ),
            recursive = TRUE,
            full.names = TRUE
          )
          
          LAT_FILES <- list.files(
            CORR_ROOT,
            pattern = paste0(
              "^FIG5B_LATENT_",
              CORR_SAMPLE_SCOPE,
              "_FEATURE_FEATURE_SPEARMAN_RHO_MATRIX\\.csv$"
            ),
            recursive = TRUE,
            full.names = TRUE
          )
          
          assert(
            length(OBS_FILES) > 0L,
            "No se encontraron matrices OBSERVED FIG5B."
          )
          
          assert(
            length(LAT_FILES) > 0L,
            "No se encontraron matrices LATENT FIG5B."
          )
          
          
          #======CORRELACION FEATURE-FEATURE SIN CONDICIONAR POR PATHWAY===========
          
          # La correlacion entre dos features es una propiedad de sus perfiles
          # observados/reconstruidos, no del pathway en el que aparecen.
          # Los multiples FIG5B se usan solo como fuentes redundantes de la misma rho.
          
          
          #======normalizar claves de features===========
          
          normalize_corr_key <- function(x) {
            
            z <- as.character(x)
            
            # por si el separador aparece literalmente como "\r" o "\n"
            z <- gsub(
              "\\r",
              "::",
              z,
              fixed = TRUE
            )
            
            z <- gsub(
              "\\n",
              "::",
              z,
              fixed = TRUE
            )
            
            # por si read.csv conserva CR/LF reales
            z <- gsub(
              "\r",
              "::",
              z,
              fixed = TRUE
            )
            
            z <- gsub(
              "\n",
              "::",
              z,
              fixed = TRUE
            )
            
            # por si R hubiera convertido el separador en punto
            z <- sub(
              "^(tx|pr|me|cl)\\.",
              "\\1::",
              z,
              ignore.case = TRUE
            )
            
            z <- gsub(
              "::+",
              "::",
              z
            )
            
            z <- trimws(
              z
            )
            
            has_prefix <- grepl(
              "^(tx|pr|me|cl)::",
              z,
              ignore.case = TRUE
            )
            
            if (any(has_prefix)) {
              
              prefix <- tolower(
                sub(
                  "^([^:]+)::.*$",
                  "\\1",
                  z[has_prefix]
                )
              )
              
              original <- sub(
                "^[^:]+::",
                "",
                z[has_prefix]
              )
              
              z[has_prefix] <- paste(
                prefix,
                original,
                sep = "::"
              )
            }
            
            z
          }
          
          
          #======normalizar claves de nuestras features finales===========
          
          CORR_FEATURES$matrix_key <- normalize_corr_key(
            paste(
              CORR_FEATURES$view_code,
              CORR_FEATURES$feature_original,
              sep = "::"
            )
          )
          
          
          #======leer una matriz FIG5B y recuperar pares disponibles===========
          
          read_rho_pairs_any <- function(
    path,
    selected_keys
          ) {
            
            d <- utils::read.csv(
              path,
              stringsAsFactors = FALSE,
              check.names = FALSE
            )
            
            #======ignorar matrices sin suficientes features===========
            
            if (ncol(d) < 3L) {
              return(NULL)
            }
            
            row_keys <- normalize_corr_key(
              d[[1L]]
            )
            
            col_keys <- normalize_corr_key(
              names(d)[-1L]
            )
            
            m <- as.matrix(
              d[
                ,
                -1L,
                drop = FALSE
              ]
            )
            
            suppressWarnings(
              storage.mode(m) <- "double"
            )
            
            
            #======evitar posibles duplicados tras normalizacion===========
            
            keep_row <- !duplicated(
              row_keys
            )
            
            keep_col <- !duplicated(
              col_keys
            )
            
            m <- m[
              keep_row,
              keep_col,
              drop = FALSE
            ]
            
            row_keys <- row_keys[
              keep_row
            ]
            
            col_keys <- col_keys[
              keep_col
            ]
            
            rownames(m) <- row_keys
            colnames(m) <- col_keys
            
            
            #======features finales presentes en esta matriz===========
            
            common <- intersect(
              selected_keys,
              intersect(
                rownames(m),
                colnames(m)
              )
            )
            
            if (length(common) < 2L) {
              return(NULL)
            }
            
            
            m <- m[
              common,
              common,
              drop = FALSE
            ]
            
            
            #======triangulo superior===========
            
            ij <- which(
              upper.tri(m),
              arr.ind = TRUE
            )
            
            if (!nrow(ij)) {
              return(NULL)
            }
            
            
            k1 <- rownames(m)[
              ij[, 1L]
            ]
            
            k2 <- colnames(m)[
              ij[, 2L]
            ]
            
            
            #======clave canonica del par===========
            
            first <- ifelse(
              k1 <= k2,
              k1,
              k2
            )
            
            second <- ifelse(
              k1 <= k2,
              k2,
              k1
            )
            
            
            out <- data.frame(
              feature1 = first,
              feature2 = second,
              
              pair_id = paste(
                first,
                second,
                sep = "\r\r"
              ),
              
              rho = as.numeric(
                m[ij]
              ),
              
              source_file = path,
              
              stringsAsFactors = FALSE
            )
            
            
            out[
              is.finite(
                out$rho
              ),
              ,
              drop = FALSE
            ]
          }
          
          
          #======leer todas las matrices disponibles===========
          
          load_corr_any <- function(
    files,
    signal_name
          ) {
            
            rows <- lapply(
              files,
              function(f) {
                
                read_rho_pairs_any(
                  path = f,
                  selected_keys = CORR_FEATURES$matrix_key
                )
              }
            )
            
            rows <- Filter(
              function(x) {
                !is.null(x) &&
                  nrow(x) > 0L
              },
              rows
            )
            
            
            assert(
              length(rows) > 0L,
              paste(
                "Ninguna matriz",
                signal_name,
                "contenia al menos dos features finales."
              )
            )
            
            
            raw <- do.call(
              rbind,
              rows
            )
            
            
            #======una misma pareja puede aparecer en muchos pathways===========
            # su correlacion deberia ser la misma; usamos mediana por robustez
            
            rho_summary <- stats::aggregate(
              rho ~
                pair_id +
                feature1 +
                feature2,
              data = raw,
              FUN = function(x) {
                
                x <- x[
                  is.finite(x)
                ]
                
                if (!length(x)) {
                  return(NA_real_)
                }
                
                stats::median(
                  x
                )
              }
            )
            
            
            n_sources <- stats::aggregate(
              source_file ~ pair_id,
              data = raw,
              FUN = function(x) {
                length(
                  unique(x)
                )
              }
            )
            
            names(
              n_sources
            )[2L] <- "n_source_matrices"
            
            
            out <- merge(
              rho_summary,
              n_sources,
              by = "pair_id",
              all.x = TRUE,
              sort = FALSE
            )
            
            out$signal <- signal_name
            
            out
          }
          
          
          #======cargar OBSERVED y LATENT===========
          
          CORR_OBS <- load_corr_any(
            OBS_FILES,
            "OBSERVED"
          )
          
          CORR_LAT <- load_corr_any(
            LAT_FILES,
            "LATENT"
          )
          
          
          cat(
            "\nPares OBSERVED recuperados:",
            nrow(CORR_OBS),
            "\n"
          )
          
          cat(
            "Pares LATENT recuperados:",
            nrow(CORR_LAT),
            "\n"
          )
          
          
          #======emparejar OBSERVED y LATENT===========
          
          CORR_PAIR_SUMMARY <- merge(
            CORR_OBS[
              ,
              c(
                "pair_id",
                "feature1",
                "feature2",
                "rho",
                "n_source_matrices"
              ),
              drop = FALSE
            ],
            
            CORR_LAT[
              ,
              c(
                "pair_id",
                "feature1",
                "feature2",
                "rho",
                "n_source_matrices"
              ),
              drop = FALSE
            ],
            
            by = c(
              "pair_id",
              "feature1",
              "feature2"
            ),
            
            all = FALSE,
            sort = FALSE,
            
            suffixes = c(
              "_observed",
              "_latent"
            )
          )
          
          
          assert(
            nrow(CORR_PAIR_SUMMARY) > 0L,
            "No hubo pares comunes entre OBSERVED y LATENT."
          )
          
          
          #======pathways compartidos como anotacion independiente===========
          
          HEAT_LINKS$pathway_pair_key <- paste(
            HEAT_LINKS$program_id,
            HEAT_LINKS$library,
            HEAT_LINKS$pathway,
            sep = "::"
          )
          
          PATHWAYS_BY_FEATURE <- split(
            HEAT_LINKS$pathway_pair_key,
            as.character(
              HEAT_LINKS$feature_id
            )
          )
          
          PATHWAYS_BY_FEATURE <- lapply(
            PATHWAYS_BY_FEATURE,
            unique
          )
          
          
          CORR_PAIR_SUMMARY$n_shared_pathways <- mapply(
            function(a, b) {
              
              pa <- PATHWAYS_BY_FEATURE[[a]]
              pb <- PATHWAYS_BY_FEATURE[[b]]
              
              if (
                is.null(pa) ||
                is.null(pb)
              ) {
                return(0L)
              }
              
              length(
                intersect(
                  pa,
                  pb
                )
              )
            },
            
            CORR_PAIR_SUMMARY$feature1,
            CORR_PAIR_SUMMARY$feature2
          )
          
          
          #======concordancia observed vs latent===========
          
          CORR_PAIR_SUMMARY$same_direction <-
            is.finite(
              CORR_PAIR_SUMMARY$rho_observed
            ) &
            is.finite(
              CORR_PAIR_SUMMARY$rho_latent
            ) &
            sign(
              CORR_PAIR_SUMMARY$rho_observed
            ) ==
            sign(
              CORR_PAIR_SUMMARY$rho_latent
            )
          
          
          #======score descriptivo conservador===========
          
          CORR_PAIR_SUMMARY$rho_concordant <-
            ifelse(
              CORR_PAIR_SUMMARY$same_direction,
              
              sign(
                CORR_PAIR_SUMMARY$rho_observed
              ) *
                pmin(
                  abs(
                    CORR_PAIR_SUMMARY$rho_observed
                  ),
                  abs(
                    CORR_PAIR_SUMMARY$rho_latent
                  )
                ),
              
              NA_real_
            )
          
          
          #======comprobar===========
          
          CORR_PAIR_SUMMARY[
            order(
              -abs(
                CORR_PAIR_SUMMARY$rho_observed
              )
            ),
            ,
            drop = FALSE
          ] |>
            head(
              20
            )
          #======labels===========
          
          LABEL_MAP <- setNames(
            CORR_FEATURES$feature_label,
            CORR_FEATURES$matrix_key
          )
          
          CORR_PAIR_SUMMARY$label1 <-
            LABEL_MAP[
              CORR_PAIR_SUMMARY$feature1
            ]
          
          CORR_PAIR_SUMMARY$label2 <-
            LABEL_MAP[
              CORR_PAIR_SUMMARY$feature2
            ]
          
          
          #======formato largo OBSERVED + LATENT===========
          
          CORR_PLOT_LONG <- rbind(
            
            data.frame(
              feature1 = CORR_PAIR_SUMMARY$feature1,
              feature2 = CORR_PAIR_SUMMARY$feature2,
              rho = CORR_PAIR_SUMMARY$rho_observed,
              signal = "OBSERVED",
              n_shared_pathways =
                CORR_PAIR_SUMMARY$n_shared_pathways,
              stringsAsFactors = FALSE
            ),
            
            data.frame(
              feature1 = CORR_PAIR_SUMMARY$feature1,
              feature2 = CORR_PAIR_SUMMARY$feature2,
              rho = CORR_PAIR_SUMMARY$rho_latent,
              signal = "LATENT / RECONSTRUCTED",
              n_shared_pathways =
                CORR_PAIR_SUMMARY$n_shared_pathways,
              stringsAsFactors = FALSE
            )
          )
          
          
          #======hacer matriz simetrica===========
          
          CORR_PLOT_LONG <- rbind(
            CORR_PLOT_LONG,
            transform(
              CORR_PLOT_LONG,
              feature1 = feature2,
              feature2 = feature1
            )
          )
          
          #======orden heatmap feature-feature===========
          
          CORR_ORDER_MODE <- "dendrogram"
          # opciones:
          # "dendrogram"
          # "omic"
          
          CORR_CLUSTER_SIGNAL <- "LATENT"
          # opciones:
          # "LATENT"
          # "OBSERVED"
          
          CORR_OMIC_ORDER <- c(
            "cl",
            "pr",
            "me",
            "tx"
          )
          
          
          #======construir matriz simetrica para clustering===========
          
          build_corr_matrix_for_order <- function(
    pair_summary,
    feature_keys,
    rho_col
          ) {
            
            m <- matrix(
              NA_real_,
              nrow = length(feature_keys),
              ncol = length(feature_keys),
              dimnames = list(
                feature_keys,
                feature_keys
              )
            )
            
            diag(m) <- 1
            
            for (
              i in seq_len(
                nrow(pair_summary)
              )
            ) {
              
              a <- pair_summary$feature1[i]
              b <- pair_summary$feature2[i]
              value <- pair_summary[[rho_col]][i]
              if (
                a %in% feature_keys &&
                b %in% feature_keys &&
                is.finite(value)
              ) {
                
                m[a, b] <- value
                m[b, a] <- value
              }
            }
            
            # solo para poder calcular el dendrograma
            # si faltase alguna correlacion
            m[
              !is.finite(m)
            ] <- 0
            
            m
          }
          
          
          #======elegir orden===========
          
          if (
            CORR_ORDER_MODE ==
            "dendrogram"
          ) {
            
            rho_col <- if (
              CORR_CLUSTER_SIGNAL ==
              "LATENT"
            ) {
              "rho_latent"
            } else {
              "rho_observed"
            }
            
            CORR_CLUSTER_MATRIX <-
              build_corr_matrix_for_order(
                pair_summary =
                  CORR_PAIR_SUMMARY,
                feature_keys =
                  CORR_FEATURES$matrix_key,
                rho_col =
                  rho_col
              )
            
            CORR_DISTANCE <- stats::as.dist(
              1 -
                CORR_CLUSTER_MATRIX
            )
            
            CORR_HCLUST <- stats::hclust(
              CORR_DISTANCE,
              method = "average"
            )
            
            FEATURE_LEVELS <-
              rownames(
                CORR_CLUSTER_MATRIX
              )[
                CORR_HCLUST$order
              ]
            
          } else {
            
            CORR_META <- CORR_FEATURES
            
            CORR_META$omic_order <- match(
              CORR_META$view_code,
              CORR_OMIC_ORDER
            )
            
            CORR_META$omic_order[
              is.na(
                CORR_META$omic_order
              )
            ] <- length(
              CORR_OMIC_ORDER
            ) + 1L
            
            CORR_META <- CORR_META[
              order(
                CORR_META$omic_order,
                -CORR_META$evidence_rank_score,
                CORR_META$feature_label
              ),
              ,
              drop = FALSE
            ]
            
            FEATURE_LEVELS <-
              CORR_META$matrix_key
          }
          
          
          #======aplicar mismo orden a ambos paneles===========
          
          CORR_PLOT_LONG$feature1 <- factor(
            CORR_PLOT_LONG$feature1,
            levels = FEATURE_LEVELS
          )
          
          CORR_PLOT_LONG$feature2 <- factor(
            CORR_PLOT_LONG$feature2,
            levels = rev(
              FEATURE_LEVELS
            )
          )
          #======texto: numero de pathways compartidos solo en OBSERVED===========
          
          CORR_PLOT_LONG$cell_label <- ifelse(
            CORR_PLOT_LONG$signal == "OBSERVED" &
              CORR_PLOT_LONG$n_shared_pathways > 0L,
            CORR_PLOT_LONG$n_shared_pathways,
            ""
          )
          
          
          #======PLOT===========
          
          PLOT_CORR_FEATURES_PILOT <-
            ggplot2::ggplot(
              CORR_PLOT_LONG,
              ggplot2::aes(
                x = feature1,
                y = feature2,
                fill = rho
              )
            ) +
            
            ggplot2::geom_tile(
              colour = "white",
              linewidth = 0.25
            ) +
            
            ggplot2::geom_text(
              ggplot2::aes(label = cell_label),
              size = HM_CELL_TEXT_SIZE
            )+
            
            ggplot2::facet_wrap(
              ~ signal,
              nrow = 1
            ) +
            
            ggplot2::scale_fill_gradient2(
              low = "#3B4CC0",
              mid = "white",
              high = "#B40426",
              midpoint = 0,
              limits = c(
                -1,
                1
              ),
              name = "Spearman rho"
            ) +
            
            ggplot2::scale_x_discrete(
              labels = LABEL_MAP
            ) +
            
            ggplot2::scale_y_discrete(
              labels = LABEL_MAP
            ) +
            
            ggplot2::labs(
              title = paste0(
                CORR_FAMILY,
                " — ",
                CORR_PATTERN,
                ": feature–feature covariance"
              ),
              subtitle = paste0(
                "Observed vs latent/reconstructed Spearman correlation.\n",
                "Numbers in OBSERVED cells = shared Q5 pathways retained in the biological story."
              ),
              x = NULL,
              y = NULL
            ) +
            
            ggplot2::theme_minimal(base_size = HM_BASE_SIZE) +
            ggplot2::theme(
              panel.grid = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_text(
                angle = HM_X_ANGLE, hjust = 1, size = HM_AXIS_X_SIZE
              ),
              axis.text.y = ggplot2::element_text(size = HM_AXIS_Y_SIZE),
              strip.text = ggplot2::element_text(face = "bold", size = HM_STRIP_TEXT_SIZE),
              legend.title = ggplot2::element_text(size = HM_LEGEND_TITLE_SIZE),
              legend.text = ggplot2::element_text(size = HM_LEGEND_TEXT_SIZE),
              legend.position = "right",
              plot.title = ggplot2::element_text(face = "bold", size = HM_TITLE_SIZE),
              plot.subtitle = ggplot2::element_text(size = 15, lineheight = 1.1,
                                                    margin = ggplot2::margin(b = 12)),
              plot.margin = ggplot2::margin(20, 25, 25, 20)
            )
          
          #======mostrar===========
          
          print(
            PLOT_CORR_FEATURES_PILOT
          )
          
          
          #======tabla para inspeccion===========
          if (isTRUE(HEATMAP_OPEN_TABLES)) View(
            CORR_PAIR_SUMMARY[
              order(
                -CORR_PAIR_SUMMARY$n_shared_pathways,
                -abs(
                  CORR_PAIR_SUMMARY$rho_concordant
                )
              ),
              ,
              drop = FALSE
            ],
            title = paste(
              HEAT_FAMILY,
              HEAT_PATTERN,
              "| Correlacion observed-latent y pathways compartidos"
            )
          )
          
          #======SIMILITUD PATHWAY-PATHWAY POR PERFIL DE FEATURES===========
          
          PATHWAY_ORDER_MODE <- "dendrogram"
          # opciones:
          # "dendrogram"
          # "library"
          
          PATHWAY_MIN_FEATURES <- 2L
          
          
          #======colapsar duplicados feature-pathway===========
          
          PATHWAY_PROFILE_LONG <- stats::aggregate(
            association_percentile ~
              feature_label +
              pathway_display +
              library,
            data = HEAT_LINKS,
            FUN = max,
            na.rm = TRUE
          )
          
          
          #======numero de features conectadas a cada pathway===========
          
          PATHWAY_SUPPORT <- stats::aggregate(
            feature_label ~
              pathway_display +
              library,
            data = PATHWAY_PROFILE_LONG,
            FUN = function(x) {
              length(
                unique(x)
              )
            }
          )
          
          names(
            PATHWAY_SUPPORT
          )[
            names(
              PATHWAY_SUPPORT
            ) == "feature_label"
          ] <- "n_linked_features"
          
          
          #======quedarnos con pathways informativos===========
          
          KEEP_PATHWAYS <- PATHWAY_SUPPORT$pathway_display[
            PATHWAY_SUPPORT$n_linked_features >=
              PATHWAY_MIN_FEATURES
          ]
          
          PATHWAY_PROFILE_LONG <-
            PATHWAY_PROFILE_LONG[
              PATHWAY_PROFILE_LONG$pathway_display %in%
                KEEP_PATHWAYS,
              ,
              drop = FALSE
            ]
          
          assert(
            length(
              unique(
                PATHWAY_PROFILE_LONG$pathway_display
              )
            ) >= 2L,
            "No hay suficientes pathways para correlacion pathway-pathway."
          )
          
          
          #======orden de features usado para construir perfiles===========
          
          FEATURES_PP <- unique(
            HEAT_SELECTED_PROGRAM[
              order(
                -HEAT_SELECTED_PROGRAM$evidence_rank_score
              ),
              "feature_label"
            ]
          )
          
          
          #======lista de pathways===========
          
          PATHWAYS_PP <- unique(
            as.character(
              PATHWAY_PROFILE_LONG$pathway_display
            )
          )
          
          
          #======matriz FEATURE x PATHWAY===========
          
          FEATURE_PATHWAY_MATRIX_PP <- matrix(
            0,
            nrow = length(
              FEATURES_PP
            ),
            ncol = length(
              PATHWAYS_PP
            ),
            dimnames = list(
              FEATURES_PP,
              PATHWAYS_PP
            )
          )
          
          
          #======rellenar asociaciones existentes===========
          
          for (
            i in seq_len(
              nrow(
                PATHWAY_PROFILE_LONG
              )
            )
          ) {
            
            f <- as.character(
              PATHWAY_PROFILE_LONG$feature_label[i]
            )
            
            p <- as.character(
              PATHWAY_PROFILE_LONG$pathway_display[i]
            )
            
            a <- as.numeric(
              PATHWAY_PROFILE_LONG$association_percentile[i]
            )
            
            if (
              f %in% rownames(
                FEATURE_PATHWAY_MATRIX_PP
              ) &&
              p %in% colnames(
                FEATURE_PATHWAY_MATRIX_PP
              ) &&
              is.finite(a)
            ) {
              
              FEATURE_PATHWAY_MATRIX_PP[
                f,
                p
              ] <- a
            }
          }
          
          
          #======eliminar pathways sin variabilidad===========
          
          PATHWAY_VARIABLE <- apply(
            FEATURE_PATHWAY_MATRIX_PP,
            2,
            function(x) {
              length(
                unique(x)
              ) >= 2L
            }
          )
          
          FEATURE_PATHWAY_MATRIX_PP <-
            FEATURE_PATHWAY_MATRIX_PP[
              ,
              PATHWAY_VARIABLE,
              drop = FALSE
            ]
          
          assert(
            ncol(
              FEATURE_PATHWAY_MATRIX_PP
            ) >= 2L,
            "Tras quitar pathways constantes quedaron menos de dos pathways."
          )
          
          
          #======correlacion PATHWAY x PATHWAY===========
          
          PATHWAY_PATHWAY_COR <- stats::cor(
            FEATURE_PATHWAY_MATRIX_PP,
            method = "spearman",
            use = "pairwise.complete.obs"
          )
          
          
          #======asegurar diagonal===========
          
          diag(
            PATHWAY_PATHWAY_COR
          ) <- 1
          
          
          #======orden pathway por dendrograma o biblioteca===========
          
          if (
            PATHWAY_ORDER_MODE ==
            "dendrogram"
          ) {
            
            PATHWAY_DISTANCE <- stats::as.dist(
              1 -
                PATHWAY_PATHWAY_COR
            )
            
            PATHWAY_HCLUST <- stats::hclust(
              PATHWAY_DISTANCE,
              method = "average"
            )
            
            PATHWAY_ORDER_PP <-
              rownames(
                PATHWAY_PATHWAY_COR
              )[
                PATHWAY_HCLUST$order
              ]
            
          } else {
            
            PATHWAY_META_PP <- unique(
              PATHWAY_PROFILE_LONG[
                ,
                c(
                  "pathway_display",
                  "library"
                ),
                drop = FALSE
              ]
            )
            
            PATHWAY_LIBRARY_ORDER_PP <- c(
              "HALLMARK",
              "REACTOME",
              "GO_BP",
              "GO_MF",
              "GO_CC"
            )
            
            PATHWAY_META_PP$library_order <- match(
              PATHWAY_META_PP$library,
              PATHWAY_LIBRARY_ORDER_PP
            )
            
            PATHWAY_META_PP$library_order[
              is.na(
                PATHWAY_META_PP$library_order
              )
            ] <- length(
              PATHWAY_LIBRARY_ORDER_PP
            ) + 1L
            
            PATHWAY_META_PP <- PATHWAY_META_PP[
              order(
                PATHWAY_META_PP$library_order,
                PATHWAY_META_PP$pathway_display
              ),
              ,
              drop = FALSE
            ]
            
            PATHWAY_ORDER_PP <- intersect(
              PATHWAY_META_PP$pathway_display,
              rownames(
                PATHWAY_PATHWAY_COR
              )
            )
          }
          
          
          #======reordenar matriz===========
          
          PATHWAY_PATHWAY_COR <-
            PATHWAY_PATHWAY_COR[
              PATHWAY_ORDER_PP,
              PATHWAY_ORDER_PP,
              drop = FALSE
            ]
          
          
          #======pasar a formato largo===========
          
          PATHWAY_PATHWAY_PLOT <- as.data.frame(
            as.table(
              PATHWAY_PATHWAY_COR
            ),
            stringsAsFactors = FALSE
          )
          
          names(
            PATHWAY_PATHWAY_PLOT
          ) <- c(
            "pathway1",
            "pathway2",
            "rho_pathway"
          )
          
          
          #======factores con el mismo orden===========
          
          PATHWAY_PATHWAY_PLOT$pathway1 <- factor(
            PATHWAY_PATHWAY_PLOT$pathway1,
            levels = PATHWAY_ORDER_PP
          )
          
          PATHWAY_PATHWAY_PLOT$pathway2 <- factor(
            PATHWAY_PATHWAY_PLOT$pathway2,
            levels = rev(
              PATHWAY_ORDER_PP
            )
          )
          
          
          #======PLOT PATHWAY x PATHWAY===========
          
          PLOT_PATHWAY_PATHWAY_SIMILARITY <-
            ggplot2::ggplot(
              PATHWAY_PATHWAY_PLOT,
              ggplot2::aes(
                x = pathway1,
                y = pathway2,
                fill = rho_pathway
              )
            ) +
            
            ggplot2::geom_tile(
              colour = "white",
              linewidth = 0.25
            ) +
            
            ggplot2::scale_fill_gradient2(
              low = "#3B4CC0",
              mid = "white",
              high = "#B40426",
              midpoint = 0,
              limits = c(
                -1,
                1
              ),
              name = "Spearman\npathway profile"
            ) +
            
            ggplot2::labs(
              title = paste0(
                HEAT_FAMILY,
                " — ",
                HEAT_PATTERN,
                ": pathway–pathway similarity"
              ),
              subtitle = paste0(
                "Each pathway is represented by its association-percentile profile\n",
                "across retained features. Spearman rho measures similarity between pathway profiles."
              ),
              x = NULL,
              y = NULL
            ) +
            ggplot2::theme_minimal(base_size = HM_BASE_SIZE) +
            ggplot2::theme(
              panel.grid = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_text(
                angle = HM_X_ANGLE, hjust = 1, vjust = 1, size = HM_AXIS_X_SIZE
              ),
              axis.text.y = ggplot2::element_text(size = HM_AXIS_Y_SIZE),
              plot.title = ggplot2::element_text(face = "bold", size = HM_TITLE_SIZE),
              plot.subtitle = ggplot2::element_text(size = 15, lineheight = 1.1,
                                                    margin = ggplot2::margin(b = 12)),
              legend.title = ggplot2::element_text(size = HM_LEGEND_TITLE_SIZE),
              legend.text = ggplot2::element_text(size = HM_LEGEND_TEXT_SIZE),
              legend.position = "right",
              plot.margin = ggplot2::margin(20, 25, 25, 20)
            )
          
          #======mostrar===========
          
          print(
            PLOT_PATHWAY_PATHWAY_SIMILARITY
          )
          
          
          #======pares de pathways mas similares===========
          
          PATHWAY_PATHWAY_PAIRS <-
            PATHWAY_PATHWAY_PLOT[
              as.character(
                PATHWAY_PATHWAY_PLOT$pathway1
              ) <
                as.character(
                  PATHWAY_PATHWAY_PLOT$pathway2
                ),
              ,
              drop = FALSE
            ]
          
          PATHWAY_PATHWAY_PAIRS <-
            PATHWAY_PATHWAY_PAIRS[
              order(
                -PATHWAY_PATHWAY_PAIRS$rho_pathway
              ),
              ,
              drop = FALSE
            ]
          
          rownames(
            PATHWAY_PATHWAY_PAIRS
          ) <- NULL
          
          #======tabla pathway-pathway opcional===========
          
          if (isTRUE(HEATMAP_OPEN_TABLES)) {
            View(
              PATHWAY_PATHWAY_PAIRS,
              title = paste(
                HEAT_FAMILY,
                HEAT_PATTERN,
                "| Pathway-pathway similarity"
              )
            )
          }
          
          #======resultados del escenario en memoria===========
          
          list(
            feature_pathway = PLOT_HEATMAP_PILOT,
            feature_covariance = PLOT_CORR_FEATURES_PILOT,
            pathway_similarity = PLOT_PATHWAY_PATHWAY_SIMILARITY,
            
            heat_links = HEAT_LINKS,
            correlation_pairs = CORR_PAIR_SUMMARY,
            pathway_pairs = PATHWAY_PATHWAY_PAIRS
          )
          
        }),
    
    error = function(e) {
      message(
        "SKIP: ",
        fam_now, " | ", pattern_now,
        " -> ", conditionMessage(e)
      )
      NULL
    }
      )
    )
  }
}

#======RESUMEN DE HEATMAPS GENERADOS===========

print(
  sapply(
    HEATMAP_RESULTS,
    function(x) {
      sum(vapply(x, function(z) !is.null(z), logical(1)))
    }
  )
)

#======PILOTO MODULOS DE FEATURES Y ANOTACIONES===========

PILOT_FAMILY <- "P1"
PILOT_PATTERN <- "SO_like_high"

PILOT_N_FEATURE_MODULES <- 3L
PILOT_N_PATHWAY_MODULES <- 3L
PILOT_CLUSTER_METHOD <- "average"

PILOT_MODULES_STATUS <- "OK"
PILOT_MODULES_SO <- NULL
PLOT_PILOT_MODULES_SO <- NULL
tryCatch({
  #======recuperar resultados ya calculados===========
  
  PILOT_SOURCE <- HEATMAP_RESULTS[[PILOT_FAMILY]][[PILOT_PATTERN]]
  
  assert(
    !is.null(PILOT_SOURCE),
    "No existen heatmaps para el piloto seleccionado."
  )
  
  PILOT_ASSOC_SOURCE <- PILOT_SOURCE$heat_links
  PILOT_CORR_SOURCE <- PILOT_SOURCE$correlation_pairs
  
  assert(
    nrow(PILOT_ASSOC_SOURCE) > 0L &&
      nrow(PILOT_CORR_SOURCE) > 0L,
    "No existen suficientes datos para construir los modulos."
  )
  
  #======adaptar datos de asociacion a las columnas reales===========
  
  ASSOC_LONG <- data.frame(
    family = PILOT_FAMILY,
    pattern = PILOT_PATTERN,
    
    feature = as.character(
      PILOT_ASSOC_SOURCE$feature_id
    ),
    
    pathway = paste(
      as.character(PILOT_ASSOC_SOURCE$library),
      as.character(PILOT_ASSOC_SOURCE$pathway),
      sep = " | "
    ),
    
    association_percentile = as.numeric(
      PILOT_ASSOC_SOURCE$association_percentile
    ),
    
    stringsAsFactors = FALSE
  )
  
  ASSOC_LONG <- stats::aggregate(
    association_percentile ~ family + pattern + feature + pathway,
    data = ASSOC_LONG,
    FUN = max
  )
  
  #======adaptar correlaciones observadas a las columnas reales===========
  
  CORR_OBS_MODULES <- data.frame(
    family = PILOT_FAMILY,
    pattern = PILOT_PATTERN,
    
    feature_1 = as.character(
      PILOT_CORR_SOURCE$feature1
    ),
    
    feature_2 = as.character(
      PILOT_CORR_SOURCE$feature2
    ),
    
    rho = as.numeric(
      PILOT_CORR_SOURCE$rho_observed
    ),
    
    stringsAsFactors = FALSE
  )
  
  CORR_OBS_MODULES <- CORR_OBS_MODULES[
    is.finite(CORR_OBS_MODULES$rho),
    ,
    drop = FALSE
  ]
  
  #======construir modulos===========
  
  PILOT_MODULES_SO <- build_feature_function_modules(
    assoc_long = ASSOC_LONG,
    corr_long = CORR_OBS_MODULES,
    
    family = PILOT_FAMILY,
    pattern = PILOT_PATTERN,
    
    n_feature_modules = PILOT_N_FEATURE_MODULES,
    n_pathway_modules = PILOT_N_PATHWAY_MODULES,
    
    feature_col = "feature",
    pathway_col = "pathway",
    assoc_col = "association_percentile",
    
    corr_feature1_col = "feature_1",
    corr_feature2_col = "feature_2",
    corr_value_col = "rho"
  )
  
  #======NO repetir aqui el ultimo heatmap del bucle===========
  
  #======HEATMAP PILOTO DE MODULOS (NO CAMBIA CLUSTERING)===========
  # Se conserva como figura adicional; no se imprime masivamente en RStudio.
  PLOT_PILOT_MODULES_SO <- plot_feature_function_module_heatmap(
    module_obj = PILOT_MODULES_SO,
    family = PILOT_FAMILY,
    pattern = PILOT_PATTERN
  )
  
  #======inspeccionar miembros de cada modulo===========
  if (!is.null(PILOT_MODULES_SO)) {
    print(PILOT_MODULES_SO$feat_module_df)
    print(PILOT_MODULES_SO$path_module_df)
  }
  # print(PILOT_MODULES_SO$module_heat)
  
}, error = function(e) {
  PILOT_MODULES_STATUS <<- paste0("ERROR: ", conditionMessage(e))
  PILOT_MODULES_SO <<- NULL
  PLOT_PILOT_MODULES_SO <<- NULL
  message("Modulo piloto omitido sin detener el resto de figuras: ",
          conditionMessage(e))
})

#======DENSIDAD DE TODAS LAS REDES (CAMBIAR SOLO ESTA LINEA)===========
NETWORK_DENSITY <- "media"   # "alta" | "media" | "baja"

# Umbrales VISUALES: no alteran features finales, MOFA ni Q5.
NETWORK_DENSITY_PRESETS <- list(
  alta  = list(ff_rho = 0.45, pp_rho = 0.60, fp_assoc = 0.15,
               ff_top_k = 4L, pp_top_k = 3L, fp_top_k = 3L),
  media = list(ff_rho = 0.60, pp_rho = 0.70, fp_assoc = 0.25,
               ff_top_k = 3L, pp_top_k = 2L, fp_top_k = 2L),
  baja  = list(ff_rho = 0.75, pp_rho = 0.85, fp_assoc = 0.45,
               ff_top_k = 2L, pp_top_k = 1L, fp_top_k = 1L)
)
stopifnot(NETWORK_DENSITY %in% names(NETWORK_DENSITY_PRESETS))
NETWORK_CFG <- NETWORK_DENSITY_PRESETS[[NETWORK_DENSITY]]

# Se trabaja con los 3 patrones principales, por pregunta (P1/P2/P3).
FRAMEWORK_PATTERNS <- c("FD_like_high", "NP_like_high", "SO_like_high")

#======FUNCIONES AUXILIARES===========
net3_topk_directed <- function(z, k) {
  if (!nrow(z)) return(z)
  sel <- unlist(lapply(split(seq_len(nrow(z)), z$from), function(i) {
    head(i[order(-z$strength[i], z$to[i])], k)
  }), use.names = FALSE)
  z[sort(unique(sel)), , drop = FALSE]
}

net3_topk_undirected <- function(z, k) {
  if (!nrow(z)) return(z)
  sel <- unlist(lapply(unique(c(z$from, z$to)), function(id) {
    i <- which(z$from == id | z$to == id)
    head(i[order(-z$strength[i], z$from[i], z$to[i])], k)
  }), use.names = FALSE)
  z[sort(unique(sel)), , drop = FALSE]
}

net3_short_path <- function(library, pathway) {
  x <- paste(library, pathway, sep = " | ")
  x <- sub("^[^|]+\\|\\s*", "", x)
  x <- sub("^(GOBP|GO_BP|GOCC|GO_CC|GOMF|GO_MF|HALLMARK|REACTOME)_", "", x)
  x <- gsub("_+", " ", x)
  vapply(x, function(z) paste(strwrap(z, width = NET_PATH_LABEL_WRAP), collapse = "\n"), character(1))
}

net3_signed_clusters <- function(ids, pairs, rho_col, k, prefix) {
  m <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
  diag(m) <- 1
  if (nrow(pairs)) {
    for (i in seq_len(nrow(pairs))) {
      a <- as.character(pairs$from[i]); b <- as.character(pairs$to[i])
      r <- suppressWarnings(as.numeric(pairs[[rho_col]][i]))
      if (a %in% ids && b %in% ids && is.finite(r)) {
        m[a, b] <- r
        m[b, a] <- r
      }
    }
  }
  n_missing <- sum(!is.finite(m[upper.tri(m)]))
  if (n_missing) {
    message(prefix, ": ", n_missing,
            " pares sin rho; distancia neutra (=1) SOLO para dendrograma.")
  }
  m[!is.finite(m)] <- 0  # ausencia: no se crean aristas a partir de estos ceros
  diag(m) <- 1
  if (length(ids) < 2L) {
    groups <- setNames(rep(1L, length(ids)), ids)
  } else {
    hc <- stats::hclust(stats::as.dist(1 - m), method = "average")
    groups <- stats::cutree(hc, k = min(as.integer(k), length(ids)))
  }
  setNames(paste0(prefix, unname(groups)), names(groups))
}

net3_profile_clusters <- function(features, pathways, raw, k) {
  # Matriz de anotacion Q5: ceros = falta de enlace Q5, no rho observado.
  mat <- matrix(0, nrow = length(features), ncol = length(pathways),
                dimnames = list(features, pathways))
  if (nrow(raw)) {
    for (i in seq_len(nrow(raw))) {
      a <- paste0("F:", raw$feature_id[i])
      b <- paste0("P:", raw$library[i], "|", raw$pathway[i])
      v <- suppressWarnings(as.numeric(raw$association_percentile[i]))
      if (a %in% features && b %in% pathways && is.finite(v)) {
        mat[a, b] <- max(mat[a, b], v)
      }
    }
  }
  if (length(pathways) == 1L) return(setNames("PM1", pathways))
  sim <- suppressWarnings(stats::cor(mat, method = "spearman"))
  if (anyNA(sim)) {
    message("Algunos perfiles Q5 son constantes: similitud 0 SOLO para clustering PM.")
    sim[!is.finite(sim)] <- 0
  }
  diag(sim) <- 1
  hc <- stats::hclust(stats::as.dist(1 - sim), method = "average")
  groups <- stats::cutree(hc, k = min(as.integer(k), length(pathways)))
  setNames(paste0("PM", unname(groups)), names(groups))
}

net3_rescale <- function(z, lim) {
  if (!length(z)) return(z)
  if (length(unique(z)) == 1L) return(rep(mean(lim), length(z)))
  lim[1L] + (z - min(z)) / diff(range(z)) * diff(lim)
}

net3_layout_side <- function(ids, edges, xlim, ylim = c(-6.3, 6.3)) {
  if (!length(ids)) stop("Se intento dibujar una capa sin nodos.")
  vertices <- data.frame(name = ids)
  graph <- igraph::graph_from_data_frame(
    edges[, c("from", "to", "strength"), drop = FALSE],
    directed = FALSE, vertices = vertices
  )
  if (length(ids) == 1L) {
    xy <- matrix(c(mean(xlim), 0), nrow = 1L)
  } else {
    xy <- igraph::layout_with_fr(
      graph, niter = 1400,
      weights = if (nrow(edges)) igraph::E(graph)$strength else NULL
    )
    xy[, 1L] <- net3_rescale(xy[, 1L], xlim)
    xy[, 2L] <- net3_rescale(xy[, 2L], ylim)
  }
  data.frame(id = igraph::V(graph)$name, x = xy[, 1L], y = xy[, 2L])
}

net3_edge_xy <- function(z, pos) {
  z$x <- pos$x[match(z$from, pos$id)]
  z$y <- pos$y[match(z$from, pos$id)]
  z$xend <- pos$x[match(z$to, pos$id)]
  z$yend <- pos$y[match(z$to, pos$id)]
  z$width_plot <- 0.30 + 1.90 * z$strength
  z
}


#======funcion red acoplada para cualquier grupo y senal===========

build_coupled_network_plot <- function(
    net_family,
    net_pattern,
    net_ff_mode = c("observed", "latent"),
    net_fm_source = "active",
    net_ff_min_rho = 0.60,
    net_pp_min_rho = 0.70,
    net_fp_min_assoc = 0.25,
    net_ff_top_k = 3L,
    net_pp_top_k = 2L,
    net_fp_top_k = 2L,
    net_require_q5_link = TRUE,
    net_k_feature_modules = 3L,
    net_k_pathway_modules = 3L,
    net_seed = 2026L
) {
  
  net_ff_mode <- match.arg(net_ff_mode)
  
  stopifnot(
    requireNamespace("igraph", quietly = TRUE),
    requireNamespace("ggplot2", quietly = TRUE),
    requireNamespace("ggrepel", quietly = TRUE),
    exists("HEATMAP_RESULTS"),
    exists("FEATURES_FINAL"),
    net_fm_source %in% c("active", "observed")
  )
  
  NET_SRC <- HEATMAP_RESULTS[[net_family]][[net_pattern]]
  if (is.null(NET_SRC)) stop("No hay HEATMAP_RESULTS para la combinacion elegida.")
  stopifnot(all(c("heat_links", "correlation_pairs", "pathway_pairs") %in% names(NET_SRC)))
  
  NET_RAW <- NET_SRC$heat_links
  NET_COR <- NET_SRC$correlation_pairs
  NET_PPAIRS <- NET_SRC$pathway_pairs
  stopifnot(nrow(NET_RAW) > 0L, nrow(NET_COR) > 0L)
  
  #======nodos de features y pathways===========
  NET_F_ALL <- FEATURES_FINAL[[net_family]]
  NET_F_ALL <- NET_F_ALL[
    NET_F_ALL$framework_pattern == net_pattern,
    c("feature_id", "feature_label", "view_code", "evidence_rank_score"),
    drop = FALSE
  ]
  NET_F_ALL <- NET_F_ALL[!duplicated(NET_F_ALL$feature_id), , drop = FALSE]
  NET_F_ALL$id <- paste0("F:", NET_F_ALL$feature_id)
  
  NET_P_ALL <- unique(NET_RAW[, c("library", "pathway", "pathway_display"), drop = FALSE])
  NET_P_ALL$id <- paste0("P:", NET_P_ALL$library, "|", NET_P_ALL$pathway)
  NET_P_ALL <- NET_P_ALL[!duplicated(NET_P_ALL$id), , drop = FALSE]
  NET_P_ALL$label <- net3_short_path(NET_P_ALL$library, NET_P_ALL$pathway)
  
  if (nrow(NET_F_ALL) < 2L || nrow(NET_P_ALL) < 2L) {
    stop("Muy pocos nodos para construir la red.")
  }
  
  #======feature-pathway===========
  NET_FP <- data.frame(
    from = paste0("F:", NET_RAW$feature_id),
    to = paste0("P:", NET_RAW$library, "|", NET_RAW$pathway),
    strength = as.numeric(NET_RAW$association_percentile),
    stringsAsFactors = FALSE
  )
  NET_FP <- NET_FP[
    is.finite(NET_FP$strength) &
      NET_FP$strength >= net_fp_min_assoc &
      NET_FP$from %in% NET_F_ALL$id &
      NET_FP$to %in% NET_P_ALL$id,
    , drop = FALSE
  ]
  if (!nrow(NET_FP)) stop("Sin enlaces Q5 tras el filtro. Baja net_fp_min_assoc.")
  NET_FP <- stats::aggregate(strength ~ from + to, NET_FP, max)
  NET_FP <- net3_topk_directed(NET_FP, net_fp_top_k)
  
  NET_F <- NET_F_ALL
  NET_P <- NET_P_ALL
  if (isTRUE(net_require_q5_link)) {
    NET_F <- NET_F[NET_F$id %in% NET_FP$from, , drop = FALSE]
    NET_P <- NET_P[NET_P$id %in% NET_FP$to, , drop = FALSE]
  }
  if (nrow(NET_F) < 2L || nrow(NET_P) < 2L) {
    stop("Muy pocos nodos tras el filtro Q5.")
  }
  
  #======correlaciones ff===========
  NET_COR$rho_observed <- as.numeric(NET_COR$rho_observed)
  NET_COR$rho_latent <- as.numeric(NET_COR$rho_latent)
  NET_COR$rho_concordant <- ifelse(
    is.finite(NET_COR$rho_observed) & is.finite(NET_COR$rho_latent) &
      sign(NET_COR$rho_observed) == sign(NET_COR$rho_latent),
    sign(NET_COR$rho_observed) * pmin(abs(NET_COR$rho_observed), abs(NET_COR$rho_latent)),
    NA_real_
  )
  NET_COR$rho_active <- switch(
    net_ff_mode,
    observed = NET_COR$rho_observed,
    latent = NET_COR$rho_latent
  )
  NET_COR$rho_cluster <- if (net_fm_source == "observed") {
    NET_COR$rho_observed
  } else {
    NET_COR$rho_active
  }
  
  NET_FF_ALL <- data.frame(
    from = paste0("F:", NET_COR$feature1),
    to = paste0("F:", NET_COR$feature2),
    rho_active = NET_COR$rho_active,
    rho_cluster = NET_COR$rho_cluster,
    stringsAsFactors = FALSE
  )
  NET_FF_ALL <- NET_FF_ALL[
    NET_FF_ALL$from %in% NET_F$id &
      NET_FF_ALL$to %in% NET_F$id &
      NET_FF_ALL$from != NET_FF_ALL$to,
    , drop = FALSE
  ]
  
  NET_FF <- NET_FF_ALL[
    is.finite(NET_FF_ALL$rho_active) &
      abs(NET_FF_ALL$rho_active) >= net_ff_min_rho,
    c("from", "to", "rho_active"),
    drop = FALSE
  ]
  NET_FF$rho <- NET_FF$rho_active
  NET_FF$strength <- abs(NET_FF$rho)
  NET_FF <- net3_topk_undirected(NET_FF, net_ff_top_k)
  
  #======pathway-pathway===========
  NET_P_DISPLAY <- as.character(NET_P$pathway_display)
  ambiguous <- duplicated(NET_P_DISPLAY) | duplicated(NET_P_DISPLAY, fromLast = TRUE)
  NET_PATH_MAP <- setNames(NET_P$id[!ambiguous], NET_P_DISPLAY[!ambiguous])
  
  NET_PP <- data.frame(
    from = unname(NET_PATH_MAP[as.character(NET_PPAIRS$pathway1)]),
    to = unname(NET_PATH_MAP[as.character(NET_PPAIRS$pathway2)]),
    rho = as.numeric(NET_PPAIRS$rho_pathway),
    stringsAsFactors = FALSE
  )
  NET_PP <- NET_PP[
    !is.na(NET_PP$from) & !is.na(NET_PP$to) &
      NET_PP$from != NET_PP$to & is.finite(NET_PP$rho) &
      abs(NET_PP$rho) >= net_pp_min_rho,
    , drop = FALSE
  ]
  NET_PP$strength <- abs(NET_PP$rho)
  NET_PP <- net3_topk_undirected(NET_PP, net_pp_top_k)
  
  #======eliminar nodos aislados===========
  NET_ACTIVE_IDS <- unique(c(
    NET_FP$from, NET_FP$to,
    NET_FF$from, NET_FF$to,
    NET_PP$from, NET_PP$to
  ))
  
  NET_F <- NET_F[NET_F$id %in% NET_ACTIVE_IDS, , drop = FALSE]
  NET_P <- NET_P[NET_P$id %in% NET_ACTIVE_IDS, , drop = FALSE]
  NET_FP <- NET_FP[NET_FP$from %in% NET_F$id & NET_FP$to %in% NET_P$id, , drop = FALSE]
  NET_FF <- NET_FF[NET_FF$from %in% NET_F$id & NET_FF$to %in% NET_F$id, , drop = FALSE]
  NET_PP <- NET_PP[NET_PP$from %in% NET_P$id & NET_PP$to %in% NET_P$id, , drop = FALSE]
  
  #======clusters===========
  NET_FM <- net3_signed_clusters(
    ids = NET_F$id,
    pairs = NET_FF_ALL,
    rho_col = "rho_cluster",
    k = net_k_feature_modules,
    prefix = "FM"
  )
  NET_F$feature_cluster <- unname(NET_FM[NET_F$id])
  
  NET_PM <- net3_profile_clusters(
    features = NET_F$id,
    pathways = NET_P$id,
    raw = NET_RAW,
    k = net_k_pathway_modules
  )
  NET_P$pathway_cluster <- unname(NET_PM[NET_P$id])
  
  #======aristas===========
  NET_FP$feature_cluster <- NET_F$feature_cluster[match(NET_FP$from, NET_F$id)]
  NET_FP$edge_class <- paste0("Q5: ", NET_FP$feature_cluster)
  NET_FF$edge_class <- ifelse(NET_FF$rho >= 0, "Features: rho +", "Features: rho -")
  NET_PP$pm1 <- NET_P$pathway_cluster[match(NET_PP$from, NET_P$id)]
  NET_PP$pm2 <- NET_P$pathway_cluster[match(NET_PP$to, NET_P$id)]
  NET_PP$edge_class <- ifelse(
    NET_PP$rho < 0,
    "Pathways: rho -",
    ifelse(NET_PP$pm1 == NET_PP$pm2,
           paste0("Pathways: ", NET_PP$pm1),
           "Pathways: entre PM")
  )
  
  #======colores===========
  NET_FM_BASE <- c(FM1 = "#00A99D", FM2 = "#E67E22", FM3 = "#D33682",
                   FM4 = "#7A3E9D", FM5 = "#1F77B4")
  NET_PM_BASE <- c(PM1 = "#F1C40F", PM2 = "#7E57C2", PM3 = "#3498DB",
                   PM4 = "#E74C3C", PM5 = "#2ECC71")
  
  NET_FM_PRESENT <- sort(unique(NET_F$feature_cluster))
  NET_PM_PRESENT <- sort(unique(NET_P$pathway_cluster))
  NET_FM_COLORS <- NET_FM_BASE[NET_FM_PRESENT]
  NET_PM_COLORS <- NET_PM_BASE[NET_PM_PRESENT]
  
  NET_EDGE_COLORS <- c(
    "Features: rho +" = "#2166AC",
    "Features: rho -" = "#B2182B",
    "Pathways: rho -" = "#D6604D",
    "Pathways: entre PM" = "#7F7F7F",
    setNames(unname(NET_FM_COLORS), paste0("Q5: ", names(NET_FM_COLORS))),
    setNames(unname(NET_PM_COLORS), paste0("Pathways: ", names(NET_PM_COLORS)))
  )
  
  NET_FILL_COLORS <- c(
    cl = "#DF444B", me = "#389848", tx = "#A45D97", pr = "#3577B8",
    NET_PM_COLORS
  )
  
  NET_F$fill_group <- NET_F$view_code
  NET_P$fill_group <- NET_P$pathway_cluster
  
  #======layout===========
  set.seed(net_seed)
  NET_POS <- rbind(
    net3_layout_side(NET_F$id, NET_FF, c(-8.3, -2.4), c(-6.2, 6.2)),
    net3_layout_side(NET_P$id, NET_PP, c( 2.4,  8.3), c(-6.2, 6.2))
  )
  
  NET_E_FP <- net3_edge_xy(NET_FP, NET_POS)
  NET_E_FF <- net3_edge_xy(NET_FF, NET_POS)
  NET_E_PP <- net3_edge_xy(NET_PP, NET_POS)
  
  NET_F <- merge(NET_F, NET_POS, by = "id", sort = FALSE)
  NET_P <- merge(NET_P, NET_POS, by = "id", sort = FALSE)
  NET_EDGE_LABELS <- c(
    "Features: rho +" = "Feature-feature rho (+)",
    "Features: rho -" = "Feature-feature rho (-)",
    "Pathways: rho -" = "Pathway-pathway rho (-)",
    "Pathways: entre PM" = "Pathway-pathway rho between PMs",
    setNames(
      paste0("Q5 links from ", names(NET_FM_COLORS), " (feature module)"),
      paste0("Q5: ", names(NET_FM_COLORS))
    ),
    setNames(
      paste0("Pathway rho within ", names(NET_PM_COLORS), " (pathway module)"),
      paste0("Pathways: ", names(NET_PM_COLORS))
    )
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_curve(
      data = NET_E_FP,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                   colour = edge_class, linewidth = width_plot),
      curvature = 0.035, alpha = 0.78, lineend = "round"
    ) +
    ggplot2::geom_curve(
      data = NET_E_FF,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                   colour = edge_class, linewidth = width_plot),
      curvature = 0.15, alpha = 0.88, lineend = "round"
    ) +
    ggplot2::geom_curve(
      data = NET_E_PP,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                   colour = edge_class, linewidth = width_plot),
      curvature = -0.17, alpha = 0.88, lineend = "round"
    ) +
    ggplot2::geom_point(
      data = NET_F,
      ggplot2::aes(x = x, y = y, fill = fill_group, size = evidence_rank_score),
      shape = 21, colour = "#1F1F1F", stroke = 0.45
    ) +
    ggplot2::geom_point(
      data = NET_P,
      ggplot2::aes(x = x, y = y, fill = fill_group),
      shape = 22, size = NET_PATHWAY_NODE_SIZE, colour = "#303039", stroke = 0.45
    )+
    ggrepel::geom_text_repel(
      data = NET_F,
      ggplot2::aes(x = x, y = y, label = feature_label),
      size = NET_FEATURE_LABEL_SIZE, box.padding = 0.45, point.padding = 0.31,
      force = 2.4, force_pull = 0.10,
      xlim = c(-10.7, -1.0),max.time = NET_GGREPEL_MAX_TIME,
      max.overlaps = Inf, seed = net_seed, min.segment.length = 0,
      segment.colour = "grey65"
    ) +
    ggrepel::geom_text_repel(
      data = NET_P,
      ggplot2::aes(x = x, y = y, label = label),
      size = NET_PATHWAY_LABEL_SIZE, box.padding = 0.47, point.padding = 0.26,
      force = 2.4, force_pull = 0.10,
      xlim = c(1.0, 10.7), max.time =NET_GGREPEL_MAX_TIME,
      max.overlaps = Inf, seed = net_seed, min.segment.length = 0,
      segment.colour = "grey65"
    ) +
    ggplot2::scale_fill_manual(
      values = NET_FILL_COLORS,
      breaks = c("cl", "me", "pr", "tx", NET_PM_PRESENT),
      labels = c(
        "CL", "ME", "PR", "TX",
        unname(NETWORK_PM_LABELS[NET_PM_PRESENT])
      ),
      name = "Nodos: omica / modulo funcional",
      drop = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = NET_EDGE_COLORS,
      breaks = names(NET_EDGE_COLORS),
      labels = NET_EDGE_LABELS[names(NET_EDGE_COLORS)],
      name = "Connections",
      drop = TRUE
    )+
    ggplot2::scale_size_continuous(
      limits = c(0, 1),
      range = NET_FEATURE_NODE_SIZE_RANGE,
      name = "Evidence rank score"
    )+
    ggplot2::scale_linewidth_identity() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::coord_cartesian(
      xlim = NET_COORD_XLIM, ylim = NET_COORD_YLIM, clip = "off"
    ) +
    ggplot2::labs(
      title = paste0(net_family, " — ", net_pattern, ": coupled network (", net_ff_mode, ")"),
      subtitle = paste0(
        "Feature correlations: ", net_ff_mode,
        ". Pathway links: Q5 profile similarity. Cross-layer: Q5 associations."
      ),
      caption = paste0(
        "Grosor: |rho| dentro de capa; percentil Q5 entre capas. ",
        "FM por correlacion ", if (net_fm_source == "observed") "observada" else net_ff_mode,
        "; PM por perfiles Q5. Correlaciones GLOBAL: descriptivas."
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_void(base_size = NET_BASE_SIZE) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(size = NET_LEGEND_TITLE_SIZE),
      legend.text = ggplot2::element_text(size = NET_LEGEND_TEXT_SIZE),
      plot.margin = grid::unit(NET_PLOT_MARGIN, "lines")
    )
  
  list(
    plot = p,
    features = NET_F,
    pathways = NET_P,
    edges_fp = NET_E_FP,
    edges_ff = NET_E_FF,
    edges_pp = NET_E_PP
  )
}

#======3) AÑADIR EL BATCH PARA TODOS LOS GRUPOS Y EN OBSERVED + LATENT===========

#======redes acopladas para todos los grupos y senales===========

NETWORK_SIGNALS <- c("observed", "latent")
NETWORK_RESULTS <- list()

for (signal_now in NETWORK_SIGNALS) {
  
  NETWORK_RESULTS[[signal_now]] <- list()
  
  for (fam_now in FAMILIES) {
    
    NETWORK_RESULTS[[signal_now]][[fam_now]] <- list()
    
    for (pattern_now in FRAMEWORK_PATTERNS) {
      
      cat("\n==============================\n")
      cat("RED ACOPLADA:", signal_now, "|", fam_now, "|", pattern_now, "\n")
      cat("==============================\n")
      
      net_try <- try(
        build_coupled_network_plot(
          net_family = fam_now,
          net_pattern = pattern_now,
          net_ff_mode = signal_now,
          net_fm_source = "active",
          net_ff_min_rho = NETWORK_CFG$ff_rho,
          net_pp_min_rho = NETWORK_CFG$pp_rho,
          net_fp_min_assoc = NETWORK_CFG$fp_assoc,
          net_ff_top_k = NETWORK_CFG$ff_top_k,
          net_pp_top_k = NETWORK_CFG$pp_top_k,
          net_fp_top_k = NETWORK_CFG$fp_top_k,
          net_require_q5_link = TRUE,
          net_k_feature_modules = NETWORK_K_FEATURE_MODULES,
          net_k_pathway_modules = NETWORK_K_PATHWAY_MODULES,
          net_seed = 2026L
        ),
        silent = TRUE
      )
      
      if (inherits(net_try, "try-error")) {
        cat("Saltado:", fam_now, "|", pattern_now, "|", signal_now, "\n")
        cat(as.character(net_try), "\n")
        NETWORK_RESULTS[[signal_now]][[fam_now]][[pattern_now]] <- NULL
      } else {
        NETWORK_RESULTS[[signal_now]][[fam_now]][[pattern_now]] <- net_try
      }
    }
  }
}

#======AUDITORIA DE LAS 18 REDES (EN MEMORIA)===========
NETWORK_AUDIT <- do.call(rbind, lapply(NETWORK_SIGNALS, function(sig) {
  do.call(rbind, lapply(FAMILIES, function(fam) {
    do.call(rbind, lapply(FRAMEWORK_PATTERNS, function(pat) {
      z <- NETWORK_RESULTS[[sig]][[fam]][[pat]]
      data.frame(signal = sig, family = fam, pattern = pat,
                 generated = !is.null(z),
                 n_features = if (is.null(z)) NA_integer_ else nrow(z$features),
                 n_pathways = if (is.null(z)) NA_integer_ else nrow(z$pathways),
                 n_ff = if (is.null(z)) NA_integer_ else nrow(z$edges_ff),
                 n_pp = if (is.null(z)) NA_integer_ else nrow(z$edges_pp),
                 n_fp = if (is.null(z)) NA_integer_ else nrow(z$edges_fp))
    }))
  }))
}))
rownames(NETWORK_AUDIT) <- NULL
print(NETWORK_AUDIT)


#======4) VISUALIZACION OPCIONAL DE REDES===========

# FALSE: no dibuja en RStudio; permite continuar con la exportacion.
# TRUE: muestra una red seleccionada en un dispositivo independiente.
SHOW_NETWORK_IN_RSTUDIO <- FALSE

NETWORK_SHOW_SIGNAL <- "observed"
NETWORK_SHOW_FAMILY <- "P1"
NETWORK_SHOW_PATTERN <- "SO_like_high"

if (isTRUE(SHOW_NETWORK_IN_RSTUDIO)) {
  
  obj <- NETWORK_RESULTS[[NETWORK_SHOW_SIGNAL]][[NETWORK_SHOW_FAMILY]][[NETWORK_SHOW_PATTERN]]
  
  if (!is.null(obj) && inherits(obj$plot, "ggplot")) {
    
    grDevices::dev.new(width = 20, height = 13)
    
    print(obj$plot)
    
  } else {
    message("La red seleccionada no esta disponible.")
  }
}


#======34 EXPORTACION FINAL: FIGURAS Y TABLAS GENERADORAS===========
# Se ejecuta al finalizar ESTE script; no requiere ejecutar el 28 por separado.
# NO refitea modelos, NO modifica la seleccion ni los parametros Q5.

#======PARAMETROS DE EXPORTACION===========
EXPORT_BASE_DIR <- file.path(SELECTION_RUN, "FINAL_STORY_FIGURES_TABLES")
EXPORT_RUN_LABEL <- paste0(
  format(Sys.time(), "%Y%m%dT%H%M%S"),
  "_density_", NETWORK_DENSITY
)
EXPORT_DPI <- 300L
EXPORT_PNG <- TRUE
EXPORT_PDF <- TRUE
EXPORT_Q5_CATALOG <- TRUE
EXPORT_MOFA_IMPORTANCE <- TRUE

#======COMPROBAR OBJETOS DISPONIBLES===========
EXPORT_REQUIRED <- c(
  "FEATURES_FINAL", "PLOTS_ALL", "PLOTS_SHARED", "PLOTS_SPECIFIC",
  "TRIPARTITE_DATA", "TRIPARTITE_PLOTS", "HEATMAP_RESULTS",
  "NETWORK_RESULTS", "NETWORK_AUDIT", "NETWORK_CFG",
  "FAMILIES", "FRAMEWORK_PATTERNS"
)
EXPORT_MISSING <- EXPORT_REQUIRED[!vapply(EXPORT_REQUIRED, exists, logical(1),
                                          envir = .GlobalEnv, inherits = TRUE)]
if (length(EXPORT_MISSING)) {
  stop("Faltan objetos antes de exportar: ", paste(EXPORT_MISSING, collapse = ", "),
       ". Ejecuta este script completo desde el inicio.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Falta ggplot2 para exportar figuras.", call. = FALSE)
}
if (!isTRUE(EXPORT_PNG) && !isTRUE(EXPORT_PDF)) {
  stop("Activa EXPORT_PNG o EXPORT_PDF.", call. = FALSE)
}

EXPORT_DIR <- file.path(EXPORT_BASE_DIR, EXPORT_RUN_LABEL)
if (dir.exists(EXPORT_DIR)) {
  i <- 2L
  while (dir.exists(paste0(EXPORT_DIR, "_", i))) i <- i + 1L
  EXPORT_DIR <- paste0(EXPORT_DIR, "_", i)
}
dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(EXPORT_DIR))

#======HELPERS PARA EXPORTAR Y REGISTRAR===========
EXPORT_FIG_LOG <- list()
EXPORT_TAB_LOG <- list()

export_table28 <- function(x, folder, filename) {
  if (is.null(x)) return(invisible(NULL))
  if (is.matrix(x)) {
    x <- data.frame(row_id = rownames(x), x, check.names = FALSE)
  }
  if (!is.data.frame(x)) {
    warning("Se omitio tabla no tabular: ", filename)
    return(invisible(NULL))
  }
  folder <- file.path(folder, "tables")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  full <- file.path(folder, filename)
  err <- tryCatch({
    utils::write.csv(x, file = full, row.names = FALSE,
                     na = "", fileEncoding = "UTF-8")
    NA_character_
  }, error = function(e) conditionMessage(e))
  EXPORT_TAB_LOG <<- append(EXPORT_TAB_LOG, list(data.frame(
    file = full,
    n_rows = nrow(x), n_columns = ncol(x),
    status = if (is.na(err)) "OK" else "ERROR",
    error = if (is.na(err)) "" else err,
    stringsAsFactors = FALSE
  )))
  if (!is.na(err)) warning("Error al exportar tabla ", full, ": ", err)
  invisible(full)
}

export_figure28 <- function(p, folder, kind, family = "", pattern = "",
                            scope = "", signal = "", width = 16, height = 10) {
  if (is.null(p) || !inherits(p, "ggplot")) {
    EXPORT_FIG_LOG <<- append(EXPORT_FIG_LOG, list(data.frame(
      kind = kind, family = family, pattern = pattern,
      scope = scope, signal = signal, status = "OMITIDA",
      png = "", pdf = "", directory = folder,
      error = "No existe el objeto ggplot para esta combinacion",
      stringsAsFactors = FALSE
    )))
    return(invisible(FALSE))
  }
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  png_path <- if (EXPORT_PNG) file.path(folder, "figure.png") else ""
  pdf_path <- if (EXPORT_PDF) file.path(folder, "figure.pdf") else ""
  errors <- character()
  if (EXPORT_PNG) {
    err <- tryCatch({
      ggplot2::ggsave(filename = png_path, plot = p,
                      width = width, height = height, units = "in",
                      dpi = EXPORT_DPI, bg = "white", limitsize = FALSE)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(err)) errors <- c(errors, paste0("PNG: ", err))
  }
  if (EXPORT_PDF) {
    err <- tryCatch({
      ggplot2::ggsave(filename = pdf_path, plot = p,
                      width = width, height = height, units = "in",
                      device = "pdf", bg = "white", limitsize = FALSE)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(err)) errors <- c(errors, paste0("PDF: ", err))
  }
  for (file_now in c(png_path, pdf_path)) {
    if (!nzchar(file_now)) next
    if (!file.exists(file_now) || is.na(file.info(file_now)$size) ||
        file.info(file_now)$size <= 0L) {
      errors <- c(errors, paste0("Archivo no creado/vacio: ", file_now))
    }
  }
  status <- if (length(errors)) "ERROR" else "OK"
  EXPORT_FIG_LOG <<- append(EXPORT_FIG_LOG, list(data.frame(
    kind = kind, family = family, pattern = pattern,
    scope = scope, signal = signal, status = status,
    png = if (EXPORT_PNG && file.exists(png_path)) png_path else "",
    pdf = if (EXPORT_PDF && file.exists(pdf_path)) pdf_path else "",
    directory = folder, error = paste(errors, collapse = " | "),
    stringsAsFactors = FALSE
  )))
  if (length(errors)) warning("Fallo una exportacion en ", folder, ": ",
                              paste(errors, collapse = " | "))
  invisible(!length(errors))
}

export_bar_height <- function(p) {
  n <- if (is.data.frame(p$data)) nrow(p$data) else 0L
  max(9, min(38, 4.5 + 0.44 * n))
}
export_trip_height <- function(z) {
  if (is.null(z)) return(12)
  nf <- if (is.data.frame(z$feature_nodes)) nrow(z$feature_nodes) else 0L
  np <- if (is.data.frame(z$pathway_nodes)) nrow(z$pathway_nodes) else 0L
  max(10, min(28, 4.5 + 0.37 * max(nf, np)))
}
export_heat_size <- function(p, kind) {
  n <- if (is.data.frame(p$data)) nrow(p$data) else 0L
  if (kind == "feature_covariance") {
    return(c(width = max(16, min(33, 12 + sqrt(n) * 0.65)),
             height = max(9, min(25, 5 + sqrt(n) * 0.6))))
  }
  if (kind == "feature_pathway") {
    nx <- if ("pathway_display" %in% names(p$data))
      length(unique(p$data$pathway_display)) else 0L
    ny <- if ("feature_label" %in% names(p$data))
      length(unique(p$data$feature_label)) else 0L
    return(c(width = max(17, min(32, 8 + 0.65 * nx)),
             height = max(11, min(32, 5 + 0.42 * ny))))
  }
  nx <- if ("pathway1" %in% names(p$data))
    length(unique(p$data$pathway1)) else 0L
  c(width = max(14, min(27, 7 + 0.60 * nx)),
    height = max(10, min(24, 4.5 + 0.60 * nx)))
}

#======00 TABLAS MAESTRAS, AUDITORIA Y PARAMETROS===========
GLOBAL_DIR <- file.path(EXPORT_DIR, "00_GLOBAL")
export_table28(FEATURES_FINAL$ALL, GLOBAL_DIR, "FINAL_FEATURES_ALL.csv")
for (fam in FAMILIES) {
  export_table28(FEATURES_FINAL[[fam]], GLOBAL_DIR,
                 paste0("FINAL_FEATURES_", fam, ".csv"))
}
export_table28(FEATURES_FINAL$SHARING, GLOBAL_DIR, "FINAL_FEATURE_SHARING.csv")
export_table28(NETWORK_AUDIT, GLOBAL_DIR, "NETWORK_AUDIT.csv")
if (exists("MASTER", inherits = TRUE))
  export_table28(MASTER, GLOBAL_DIR, "EVIDENCE_A_ALL_MODEL_VOTES.csv")
if (isTRUE(EXPORT_MOFA_IMPORTANCE) && exists("MOFA_IMPORTANCE", inherits = TRUE))
  export_table28(MOFA_IMPORTANCE, GLOBAL_DIR, "MOFA_NATIVE_IMPORTANCE.csv")
if (isTRUE(EXPORT_Q5_CATALOG) && exists("Q5_CATALOG_TRIP", inherits = TRUE))
  export_table28(Q5_CATALOG_TRIP, GLOBAL_DIR, "Q5_FROZEN_ANNOTATION_CATALOG.csv")

NETWORK_DENSITY_TABLE <- do.call(rbind, lapply(
  names(NETWORK_DENSITY_PRESETS), function(nm) {
    d <- NETWORK_DENSITY_PRESETS[[nm]]
    data.frame(density = nm, active = nm == NETWORK_DENSITY,
               ff_rho = d$ff_rho, pp_rho = d$pp_rho,
               fp_assoc = d$fp_assoc, ff_top_k = d$ff_top_k,
               pp_top_k = d$pp_top_k, fp_top_k = d$fp_top_k)
  }))
export_table28(NETWORK_DENSITY_TABLE, GLOBAL_DIR, "NETWORK_DENSITY_PRESETS.csv")

#======01 BARPLOTS: TODAS, COMPARTIDAS Y ESPECIFICAS===========
for (fam in FAMILIES) {
  for (direction in DIRECTION_LEVELS) {
    for (scope in c("ALL", "SHARED")) {
      p <- if (scope == "ALL") PLOTS_ALL[[fam]][[direction]] else
        PLOTS_SHARED[[fam]][[direction]]
      folder <- file.path(EXPORT_DIR, fam, "BARPLOTS", scope, direction)
      if (!is.null(p)) {
        export_table28(p$data, folder, "plot_data.csv")
      }
      export_figure28(p, folder, "barplot", fam, direction, scope,
                      width = 16, height = if (is.null(p)) 9 else export_bar_height(p))
    }
  }
  p <- PLOTS_SPECIFIC[[fam]]
  folder <- file.path(EXPORT_DIR, fam, "BARPLOTS", "SPECIFIC", "FD_NP_SO_TRELLIS")
  if (!is.null(p)) export_table28(p$data, folder, "plot_data.csv")
  export_figure28(p, folder, "barplot", fam, "FD_NP_SO", "SPECIFIC",
                  width = 17, height = if (is.null(p)) 9 else export_bar_height(p))
}

#======02 REDES TRIPARTITAS: ALL / SHARED / SPECIFIC===========
for (fam in FAMILIES) {
  for (pat in FRAMEWORK_PATTERNS) {
    for (scope in c("all", "shared", "specific")) {
      z <- TRIPARTITE_DATA[[fam]][[pat]][[scope]]
      ps <- TRIPARTITE_PLOTS[[toupper(scope)]]
      p <- if (is.null(ps)) NULL else ps[[fam]][[pat]]
      folder <- file.path(EXPORT_DIR, fam, pat, "TRIPARTITE", toupper(scope))
      if (!is.null(z)) {
        for (nm in c("feature_nodes", "program_nodes", "pathway_nodes",
                     "fp_edges", "pp_edges")) {
          export_table28(z[[nm]], folder, paste0(nm, ".csv"))
        }
      }
      export_figure28(p, folder, "tripartite", fam, pat, toupper(scope),
                      width = 22, height = export_trip_height(z))
    }
  }
}

#======03 HEATMAPS: FIGURA, DATOS EXACTOS Y TABLAS FUENTE===========
HEAT_KINDS <- c("feature_pathway", "feature_covariance", "pathway_similarity")
for (fam in FAMILIES) {
  for (pat in FRAMEWORK_PATTERNS) {
    z <- HEATMAP_RESULTS[[fam]][[pat]]
    src_folder <- file.path(EXPORT_DIR, fam, pat, "HEATMAPS", "00_SOURCE_TABLES")
    if (!is.null(z)) {
      export_table28(z$heat_links, src_folder, "feature_pathway_Q5_links.csv")
      export_table28(z$correlation_pairs, src_folder,
                     "feature_feature_observed_latent_rho.csv")
      export_table28(z$pathway_pairs, src_folder,
                     "pathway_pathway_Q5_profile_rho.csv")
    }
    for (kind in HEAT_KINDS) {
      p <- if (is.null(z)) NULL else z[[kind]]
      folder <- file.path(EXPORT_DIR, fam, pat, "HEATMAPS", kind)
      if (!is.null(p) && is.data.frame(p$data))
        export_table28(p$data, folder, "plot_data.csv")
      wh <- if (is.null(p)) {
        c(width = 16, height = 10)
      } else {
        export_heat_size(p, kind)
      }
      
      # Aumentar nombres de features en los heatmaps de correlación
      if (!is.null(p) && kind == "feature_covariance") {
        
        p <- p +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(
              size = 22,
              angle = 75,
              hjust = 1,
              vjust = 1,
              face = "bold"
            ),
            axis.text.y = ggplot2::element_text(
              size = 22,
              face = "bold"
            ),
            strip.text = ggplot2::element_text(
              size = 16,
              face = "bold"
            )
          )
        
        wh <- c(width = 22, height = 17)
      }
      
      export_figure28(
        p, folder, kind, fam, pat,
        width = unname(wh["width"]),
        height = unname(wh["height"])
      )
    }
  }
}

#======04 REDES ACOPLADAS: OBSERVED / LATENT===========
for (sig in c("observed", "latent")) {
  for (fam in FAMILIES) {
    for (pat in FRAMEWORK_PATTERNS) {
      z <- NETWORK_RESULTS[[sig]][[fam]][[pat]]
      folder <- file.path(EXPORT_DIR, fam, pat, "COUPLED_NETWORK", toupper(sig))
      if (!is.null(z)) {
        for (nm in c("features", "pathways", "edges_fp", "edges_ff", "edges_pp")) {
          export_table28(z[[nm]], folder, paste0(nm, ".csv"))
        }
        param <- data.frame(
          signal = sig, family = fam, framework_pattern = pat,
          network_density = NETWORK_DENSITY,
          ff_min_abs_rho = NETWORK_CFG$ff_rho,
          pp_min_abs_rho = NETWORK_CFG$pp_rho,
          fp_min_q5_percentile = NETWORK_CFG$fp_assoc,
          ff_top_k = NETWORK_CFG$ff_top_k,
          pp_top_k = NETWORK_CFG$pp_top_k,
          fp_top_k = NETWORK_CFG$fp_top_k,
          fm_clustering_signal = sig,
          n_fm_target = NETWORK_K_FEATURE_MODULES,
          n_pm_target = NETWORK_K_PATHWAY_MODULES,
          layout_seed = 2026L,
          stringsAsFactors = FALSE
        )
        export_table28(param, folder, "network_parameters.csv")
      }
      export_figure28(if (is.null(z)) NULL else z$plot, folder,
                      "coupled_network", fam, pat, signal = sig,
                      width = 24, height = 15)
    }
  }
}

#======05 FIGURAS ADICIONALES EXISTENTES: PILOTOS===========
# Preview tripartito (coincide con una combinacion ALL, pero se conserva de modo
# independiente para que ningun objeto grafico creado quede sin registrar).
if (exists("PLOT_TRIPARTITE_PREVIEW", inherits = TRUE)) {
  preview_folder <- file.path(EXPORT_DIR, "00_GLOBAL", "PILOT_PREVIEW_TRIPARTITE")
  if (exists("TRIP_DATA_PREVIEW", inherits = TRUE) &&
      !is.null(TRIP_DATA_PREVIEW)) {
    for (nm in c("feature_nodes", "program_nodes", "pathway_nodes",
                 "fp_edges", "pp_edges")) {
      export_table28(TRIP_DATA_PREVIEW[[nm]], preview_folder, paste0(nm, ".csv"))
    }
  }
  export_figure28(PLOT_TRIPARTITE_PREVIEW, preview_folder,
                  "tripartite_preview", TRIP_PREVIEW_FAMILY,
                  TRIP_PREVIEW_PATTERN, TRIP_PREVIEW_SCOPE,
                  width = 22,
                  height = if (exists("TRIP_DATA_PREVIEW", inherits = TRUE))
                    export_trip_height(TRIP_DATA_PREVIEW) else 12)
}

# Heatmap de modulos feature/funcion; ahora se genera con la firma correcta
# de plot_feature_function_module_heatmap(module_obj=...).
module_folder <- file.path(EXPORT_DIR, PILOT_FAMILY, PILOT_PATTERN,
                           "HEATMAPS", "feature_function_modules")
export_figure28(PLOT_PILOT_MODULES_SO, module_folder,
                "feature_function_modules", PILOT_FAMILY, PILOT_PATTERN,
                width = 15, height = 10)
if (!is.null(PILOT_MODULES_SO)) {
  for (nm in c("module_heat", "feat_module_df", "path_module_df",
               "assoc_cells", "assoc_mat", "feat_corr_mat",
               "path_profile_corr")) {
    export_table28(PILOT_MODULES_SO[[nm]], module_folder, paste0(nm, ".csv"))
  }
}

# Resumen del inventario por tipo: ayuda a comprobar barras, heatmaps y redes.
if (length(EXPORT_FIG_LOG)) {
  figure_inventory <- do.call(rbind, EXPORT_FIG_LOG)
  inventory <- as.data.frame(table(figure_inventory$kind,
                                   figure_inventory$status),
                             stringsAsFactors = FALSE)
  names(inventory) <- c("figure_kind", "status", "n")
  inventory <- inventory[inventory$n > 0L, , drop = FALSE]
  export_table28(inventory, file.path(EXPORT_DIR, "00_GLOBAL"),
                 "FIGURE_INVENTORY_BY_KIND.csv")
}

#======05 MANIFIESTOS Y PROVENIENCIA===========
EXPORT_FIGURE_MANIFEST <- do.call(rbind, EXPORT_FIG_LOG)
EXPORT_TABLE_MANIFEST <- do.call(rbind, EXPORT_TAB_LOG)
rownames(EXPORT_FIGURE_MANIFEST) <- NULL
rownames(EXPORT_TABLE_MANIFEST) <- NULL
export_table28(EXPORT_FIGURE_MANIFEST, GLOBAL_DIR, "FIGURE_MANIFEST.csv")
export_table28(EXPORT_TABLE_MANIFEST, GLOBAL_DIR, "TABLE_MANIFEST.csv")

README_EXPORT <- c(
  "FIGURAS Y TABLAS PARA LA HISTORIA BIOLOGICA MULTIMODAL",
  paste0("Exportado: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Selection run: ", SELECTION_RUN),
  if (exists("BIO_RUN", inherits = TRUE)) paste0("Biological run: ", BIO_RUN),
  if (exists("PANEL10_RUN", inherits = TRUE)) paste0("Predictive panel run: ", PANEL10_RUN),
  if (exists("SGL16_RUN", inherits = TRUE)) paste0("SGL run: ", SGL16_RUN),
  if (exists("PLS17_RUN", inherits = TRUE)) paste0("PLS run: ", PLS17_RUN),
  paste0("Network density: ", NETWORK_DENSITY),
  "",
  "ESTRUCTURA: 00_GLOBAL/ y carpetas P1/P2/P3/ con BARPLOTS/ y <pattern>/.",
  "Dentro de <pattern>: TRIPARTITE, HEATMAPS y COUPLED_NETWORK.",
  "Cada carpeta de figura: figure.png, figure.pdf y tables/*.csv cuando hay datos.",
  "Los plot_data.csv contienen los datos de entrada del ggplot correspondiente.",
  "Los nodos/aristas exportados son los utilizados para representar cada red.",
  "HEATMAPS/00_SOURCE_TABLES conserva los datos compartidos entre los 3 heatmaps.",
  "Las rutas Q5 representan asociaciones multimodales, no pertenencia metabolica directa.",
  "rho observado/latente: las fuentes FIG5B se etiquetaron GLOBAL aguas arriba.",
  "Por ello esta visualizacion no constituye validacion TEST independiente.",
  "Los umbrales de las redes son SOLO de visualizacion; las features finales no cambian.",
  "FM y PM son clusters descriptivos; su numeracion no implica equivalencia entre grupos.",
  "",
  "Consultas: 00_GLOBAL/tables/FIGURE_MANIFEST.csv y TABLE_MANIFEST.csv.",
  "Las combinaciones no generadas aparecen como OMITIDA en FIGURE_MANIFEST.csv.",
  "Se exporta tambien el heatmap piloto de modulos feature x funcion, cuando existen datos.",
  paste0("Estado del piloto de modulos: ", PILOT_MODULES_STATUS),
  "El preview tripartito se exporta adicionalmente, aun si coincide con una figura ALL."
)
writeLines(README_EXPORT, con = file.path(EXPORT_DIR, "README_EXPORT.txt"),
           useBytes = TRUE)
writeLines(capture.output(utils::sessionInfo()),
           file.path(EXPORT_DIR, "SESSION_INFO.txt"))

cat("\n=====================================================\n")
cat("EXPORTACION FINAL COMPLETADA\n")
cat("Carpeta: ", EXPORT_DIR, "\n", sep = "")
cat("Figuras OK: ", sum(EXPORT_FIGURE_MANIFEST$status == "OK"), "\n", sep = "")
cat("Figuras OMITIDAS: ", sum(EXPORT_FIGURE_MANIFEST$status == "OMITIDA"), "\n", sep = "")
cat("Figuras ERROR: ", sum(EXPORT_FIGURE_MANIFEST$status == "ERROR"), "\n", sep = "")
cat("Tablas OK: ", sum(EXPORT_TABLE_MANIFEST$status == "OK"), "\n", sep = "")
cat("=====================================================\n")
# Auditoria fisica: no basta con ggplot generado; PNG/PDF han de existir.
FIGURE_FILES <- unlist(EXPORT_FIGURE_MANIFEST[c("png", "pdf")], use.names = FALSE)
FIGURE_FILES <- FIGURE_FILES[!is.na(FIGURE_FILES) & nzchar(FIGURE_FILES)]
FIGURE_FILES_MISSING <- FIGURE_FILES[
  !file.exists(FIGURE_FILES) |
    is.na(file.info(FIGURE_FILES)$size) |
    file.info(FIGURE_FILES)$size <= 0
]
if (length(FIGURE_FILES_MISSING)) {
  writeLines(FIGURE_FILES_MISSING,
             file.path(GLOBAL_DIR, "MISSING_OR_EMPTY_FIGURES.txt"))
  warning("Archivos de figuras faltantes o vacios: ",
          length(FIGURE_FILES_MISSING))
}
cat("Archivos fisicos PNG/PDF: ", length(FIGURE_FILES), "\n", sep = "")
cat("Faltantes/vacios: ", length(FIGURE_FILES_MISSING), "\n", sep = "")
if (any(EXPORT_FIGURE_MANIFEST$status == "OMITIDA")) {
  message("Hay combinaciones sin datos o sin plot: revisa FIGURE_MANIFEST.csv. ",
          "No se han creado figuras ficticias.")
}
if (any(EXPORT_FIGURE_MANIFEST$status == "ERROR") ||
    any(EXPORT_TABLE_MANIFEST$status == "ERROR")) {
  warning("Hay errores de exportacion: revisa los manifiestos en 00_GLOBAL/tables/.")
}
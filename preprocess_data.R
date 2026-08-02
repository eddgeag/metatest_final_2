##====Carga librerias=====
library(readxl)
library(dplyr)
library(ggplot2)
library(zoo)
library(ggrepel)
# library(mixOmics)
# library(tidyr)

# ============================================================
# CONFIGURACIÓN PARAMETRIZABLE DEL ANÁLISIS
# Se puede controlar desde bash con variables de entorno.
# Se controla desde bash mediante run_pipeline_from_template.sh.
# ANALYSIS_NAME y CLINICAL_INCLUDE son obligatorios.
# ============================================================

get_env_chr <- function(name, default = NULL, required = FALSE) {
  val <- Sys.getenv(name, unset = NA_character_)
  
  if (is.na(val) || !nzchar(val)) {
    if (required) {
      stop("Variable de entorno obligatoria no definida: ", name)
    }
    return(default)
  }
  
  val
}

get_env_num <- function(name, default) {
  val <- Sys.getenv(name, unset = NA_character_)
  
  if (is.na(val) || !nzchar(val)) {
    return(default)
  }
  
  out <- suppressWarnings(as.numeric(val))
  
  if (!is.finite(out)) {
    stop("Variable de entorno numérica inválida: ", name, " = ", val)
  }
  
  out
}

parse_env_vec <- function(x, sep = "\\|") {
  if (is.null(x) || is.na(x) || !nzchar(x)) {
    return(character(0))
  }
  
  out <- trimws(strsplit(x, sep)[[1]])
  out <- out[nzchar(out)]
  unique(out)
}

ascii_feature_names <- function(x) {
  x0 <- as.character(x)
  
  x1 <- iconv(
    x0,
    from = "",
    to = "ASCII//TRANSLIT",
    sub = "_"
  )
  
  x1 <- trimws(x1)
  x1 <- gsub("[\r\n\t]+", " ", x1)
  x1 <- gsub(" +", " ", x1)
  
  bad <- is.na(x1) | !nzchar(x1)
  if (any(bad)) {
    x1[bad] <- paste0("feature_", which(bad))
  }
  
  make.unique(x1, sep = "__dup")
}



ANALYSIS_NAME <- get_env_chr("ANALYSIS_NAME", required = TRUE)
OUTDIR_BASE   <- get_env_chr("OUTDIR_BASE", ".")

Q_TX <- get_env_num("Q_TX", 0.97)
Q_PR <- get_env_num("Q_PR", 0.40)
Q_ME <- get_env_num("Q_ME", 0.97)
Q_CL <- get_env_num("Q_CL", 0.70)

OUTPUT_DIR <- file.path(OUTDIR_BASE, ANALYSIS_NAME)
PLOTS_DIR  <- file.path(OUTPUT_DIR, "plots")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR,  recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Variables clínicas candidatas
# OJO: por ahora conviene usar variables numéricas positivas,
# porque process_clinicos aplica log() y luego escala.
# ============================================================

CLINICAL_CANDIDATES <- c(
  "Edad",
  "GGT",
  "Glu",
  "Col total",
  "LDL",
  "HDL",
  "TG",
  "PCR",
  "Insulina",
  "HOMA-IR",
  "QUICKI",
  "TG-INDEX",
  "FLI",
  "TG/HDL",
  "Peso",
  "Altura",
  "IMC",
  "Ratio cintura_0M",
  "Ratio-cc",
  "%BF-DEXA",
  "CUN-BAE",
  "VAT",
  "VAT (%)",
  "Masa_grasa_0M",
  "%VAT_grasa_Total_0M",
  "kg Ginoide",
  "% Ginoide",
  "% Ginoide_grasa_Total",
  "kg Androide",
  "% Androide",
  "% Androide_grasa_Total",
  "Magra (%)",
  "Magra (kg)",
  "Ratio and/gin",
  "FFM- 0M (masa magra total + masa ósea total)",
  "Sis",
  "Dias"
)

# Default: reproduce tu análisis FLI actual
CLINICAL_DEFAULT <- c(
  "GGT",
  "Glu",
  "Col total",
  "LDL",
  "HDL",
  "TG",
  "PCR",
  "Insulina",
  "HOMA-IR",
  "QUICKI",
  "TG-INDEX",
  "TG/HDL",
  "FLI"
)

CLINICAL_INCLUDE_RAW <- get_env_chr(
  "CLINICAL_INCLUDE",
  required = TRUE
)

CLINICAL_EXCLUDE_RAW <- get_env_chr(
  "CLINICAL_EXCLUDE",
  default = ""
)


CLINICAL_INCLUDE <- parse_env_vec(CLINICAL_INCLUDE_RAW)
CLINICAL_EXCLUDE <- parse_env_vec(CLINICAL_EXCLUDE_RAW)

if (length(CLINICAL_INCLUDE) == 1 && toupper(CLINICAL_INCLUDE) == "ALL") {
  CLINICAL_VARS <- CLINICAL_CANDIDATES
} else {
  CLINICAL_VARS <- CLINICAL_INCLUDE
}

# La exclusión siempre gana sobre la inclusión
CLINICAL_VARS <- unique(setdiff(CLINICAL_VARS, CLINICAL_EXCLUDE))

missing_candidates <- setdiff(CLINICAL_VARS, CLINICAL_CANDIDATES)

if (length(missing_candidates) > 0) {
  stop(
    "Estas variables clínicas no están en CLINICAL_CANDIDATES: ",
    paste(missing_candidates, collapse = ", ")
  )
}

if (!length(CLINICAL_VARS)) {
  stop("CLINICAL_VARS quedó vacío. Revisa CLINICAL_INCLUDE / CLINICAL_EXCLUDE.")
}

cat("\n============================================================\n")
cat("ANÁLISIS:", ANALYSIS_NAME, "\n")
cat("OUTDIR_BASE:", OUTDIR_BASE, "\n")
cat("OUTPUT_DIR:", OUTPUT_DIR, "\n")
cat("PLOTS_DIR:", PLOTS_DIR, "\n")
cat("Q_TX:", Q_TX, "\n")
cat("Q_PR:", Q_PR, "\n")
cat("Q_ME:", Q_ME, "\n")
cat("Q_CL:", Q_CL, "\n")
cat("CLINICAL_INCLUDE_RAW:", CLINICAL_INCLUDE_RAW, "\n")
cat("CLINICAL_EXCLUDE_RAW:", CLINICAL_EXCLUDE_RAW, "\n")
cat("Variables clínicas finales usadas:\n")
print(CLINICAL_VARS)
cat("============================================================\n")

##=== Funciones====
calcular_mahalanobis_general <- function(datos,
                                         variables = NULL,
                                         umbral_chi2 = 0.975) {
  #' @param datos: Dataframe con las variables numéricas.
  #' @param variables: Vector con nombres de columnas a usar (si NULL, usa todas las numéricas).
  #' @param umbral_chi2: Percentil para el umbral de outlier (ej. 0.975 para 97.5%).
  
  # 1. Seleccionar variables
  if (is.null(variables)) {
    datos <- datos %>% select_if(is.numeric)  # Usar solo columnas numéricas
  } else {
    datos <- datos %>% select(all_of(variables))  # Usar variables especificadas
  }
  
  # 2. Verificar que hay datos válidos
  if (ncol(datos) == 0)
    stop("No hay variables numéricas para calcular Mahalanobis.")
  if (nrow(datos) < 2)
    stop("Se requieren al menos 2 observaciones.")
  
  # 3. Calcular media y matriz de covarianza (con manejo de NA/singularidades)
  media <- tryCatch(
    colMeans(datos, na.rm = TRUE),
    error = function(e)
      rep(NA, ncol(datos))
  )
  
  cov_matrix <- tryCatch(
    cov(datos, use = "complete.obs"),
    error = function(e) {
      warning("Matriz de covarianza singular. Se usará una matriz diagonal.")
      diag(ncol(datos)) * 1e6  # Matriz diagonal grande como fallback
    }
  )
  
  # 4. Calcular distancias
  dist_mahal <- tryCatch(
    mahalanobis(datos, center = media, cov = cov_matrix),
    error = function(e) {
      warning("Error al calcular distancias. Verifique los datos.")
      rep(NA, nrow(datos))
    }
  )
  
  # 5. Identificar outliers
  umbral <- qchisq(umbral_chi2, df = ncol(datos))
  outliers <- dist_mahal > umbral
  
  # 6. Resultado
  resultados <- data.frame(
    indice = rownames(datos),
    dist_mahal = dist_mahal,
    outlier = outliers,
    umbral_chi2 = umbral
  )
  
  return(resultados)
}
split_matrix <- function(mat, block = 6) {
  # índices de inicio de cada bloque
  starts <- seq(1, ncol(mat), by = block)
  # para cada inicio, extrae columnas i:(i+block-1)
  lapply(starts, function(i) {
    mat[, i:min(i + block - 1, ncol(mat)), drop = FALSE]
  })
}

fun_process_transcriptomics <- function(path) {
  raw <- suppressMessages(read_excel(path, sheet = 1, col_names = FALSE))
  
  # 1) Metadata de features
  feat_meta_names <- raw %>% slice(7) %>% dplyr::select(1:18) %>% unlist() %>% as.character()
  feat_meta <- raw %>%
    slice(8:n()) %>%
    dplyr::select(1:18) %>%
    as.data.frame(stringsAsFactors = FALSE)
  colnames(feat_meta) <- feat_meta_names
  
  # 2) Metadata de muestras
  sample_meta_names <- raw %>% slice(1:6) %>% pull(19) %>% as.character()
  sample_meta_vals <- raw %>%
    slice(1:6) %>%
    dplyr::select(20:60) %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(sample_meta_vals) <- sample_meta_names
  sample_meta_vals <- t(sample_meta_vals)
  sample_meta_vals <- as.data.frame(sample_meta_vals, stringsAsFactors = FALSE)
  
  code_row <- which(tolower(sample_meta_names) %in% c("codigo final", "código final"))
  if (length(code_row) == 0) stop("No se encontró columna 'CÓDIGO FINAL' en metadata de muestras.")
  rownames(sample_meta_vals) <- sample_meta_vals[, code_row[1]]
  col_code <- colnames(sample_meta_vals)[code_row[1]]
  
  # 3) Matriz de expresión por feature
  data_mat <- raw %>% slice(8:n()) %>% select(20:60) %>% as.data.frame()
  colnames(data_mat) <- sample_meta_vals[[col_code]]
  rownames(data_mat) <- rownames(feat_meta)
  
  # Depuración según tu flujo original
  data_mat <- bind_cols(GeneSymbol = feat_meta$GeneSymbol, data_mat)
  auxlog <- complete.cases(data_mat)
  data_mat <- data_mat[auxlog, ]
  feat_meta <- feat_meta[auxlog, ]
  data_mat$GeneSymbol <- feat_meta$PrimaryAccession
  rownames(data_mat) <- feat_meta$FeatureNum
  expr_feat <- data_mat[, -1, drop = FALSE]
  
  grupos <- as.factor(gsub("[-_]\\d+$", "", sample_meta_vals[[col_code]]))
  sample_meta_vals$grupo <- grupos
  
  # 4) Anotación y best_id
  nz <- function(x) { x <- trimws(as.character(x)); x[x == ""] <- NA; x }
  anno <- feat_meta
  anno$EnsemblID        <- nz(anno$EnsemblID)
  anno$EntrezGeneID     <- nz(anno$EntrezGeneID)
  anno$RefSeqAccession  <- nz(anno$RefSeqAccession)
  anno$GeneSymbol       <- nz(anno$GeneSymbol)
  
  anno$best_id <- dplyr::coalesce(anno$EnsemblID,
                                  anno$EntrezGeneID,
                                  anno$RefSeqAccession,
                                  anno$GeneSymbol)
  anno$id_type <- dplyr::case_when(
    !is.na(anno$EnsemblID)       ~ "EnsemblID",
    !is.na(anno$EntrezGeneID)    ~ "EntrezGeneID",
    !is.na(anno$RefSeqAccession) ~ "RefSeqAccession",
    !is.na(anno$GeneSymbol)      ~ "GeneSymbol",
    TRUE                         ~ "NA"
  )
  anno <- anno[match(rownames(expr_feat), anno$FeatureNum), ]
  
  # Excluir controles si existieran
  if ("ControlType" %in% colnames(anno)) {
    keep <- is.na(anno$ControlType) | anno$ControlType == 0
    expr_feat <- expr_feat[keep, , drop = FALSE]
    anno <- anno[keep, , drop = FALSE]
  }
  
  coverage <- list(
    by_id_type = table(anno$id_type, useNA = "ifany"),
    prop_annotated = mean(!is.na(anno$best_id))
  )
  
  # 5) Colapso por gen/transcrito usando mayor IQR
  valid <- !is.na(anno$best_id)
  expr_valid <- expr_feat[valid, , drop = FALSE]
  anno_valid <- anno[valid, , drop = FALSE]
  
  if (nrow(expr_valid) > 0) {
    iqr_feat <- apply(expr_valid, 1, IQR, na.rm = TRUE)
    idx_list <- split(seq_len(nrow(expr_valid)), anno_valid$best_id)
    pick <- vapply(idx_list, function(ix) ix[which.max(iqr_feat[ix])], integer(1))
    
    expr_gene <- expr_valid[pick, , drop = FALSE]
    rn_best <- anno_valid$best_id[pick]
    rownames(expr_gene) <- rn_best
    
    # columnas legibles
    expr_gene <- cbind(
      BestID     = rn_best,
      GeneSymbol = anno_valid$GeneSymbol[pick],
      expr_gene
    )
    
    anno_gene <- anno_valid[pick, c("best_id","id_type","GeneSymbol","GeneName",
                                    "EnsemblID","EntrezGeneID","RefSeqAccession",
                                    "GenbankAccession","Cytoband"), drop = FALSE]
  } else {
    expr_gene <- expr_valid
    anno_gene <- anno_valid[, c("best_id","id_type","GeneSymbol","GeneName",
                                "EnsemblID","EntrezGeneID","RefSeqAccession",
                                "GenbankAccession","Cytoband"), drop = FALSE]
  }
  entrez_vec <- anno_gene$EntrezGeneID
  okE <- !is.na(entrez_vec) & entrez_vec != ""
  
  dg <- expr_gene[okE, -(1:2), drop = FALSE]  # expresión (sin BestID/GeneSymbol)
  ag <- anno_gene[okE, , drop = FALSE]        # anotación alineada
  
  iqr_rows <- apply(dg, 1, IQR, na.rm = TRUE)
  idx_list <- split(seq_len(nrow(dg)), entrez_vec[okE])
  pickE <- vapply(idx_list, function(ix) ix[which.max(iqr_rows[ix])], integer(1))
  
  data_by_entrez <- dg[pickE, , drop = FALSE]
  anno_by_entrez <- ag[pickE, , drop = FALSE]
  rownames(data_by_entrez) <- anno_by_entrez$EntrezGeneID
  rownames(anno_by_entrez) <- anno_by_entrez$EntrezGeneID
  # 6) Salida
  list(
    feature_metadata = feat_meta,
    samples_metadata = sample_meta_vals,
    data_by_feature  = expr_feat,
    anno_by_feature  = anno,
    data_by_gene     = expr_gene[,-c(1:2)],
    anno_by_gene     = anno_gene,
    # >>> NUEVO <<<
    data_by_entrez   = data_by_entrez,   # filas = EntrezID únicos
    anno_by_entrez   = anno_by_entrez,   # filas = EntrezID únicos
    coverage         = coverage,
    id_priority      = c("EnsemblID","EntrezGeneID","RefSeqAccession","GeneSymbol"),
    code_field       = col_code
  )
  
}


prepare_transcriptomics_for_model <- function(
    trans_obj, metadata, y,
    fixed_keep_ids = NULL, fixed_train_ids = NULL, fixed_test_ids = NULL,
    p_train = 0.7, seed = 123, min_test_per_class = 1,
    validate_with_bitr = FALSE, orgdb = "org.Hs.eg.db"   # <- opcional
){
  # --- 0) Construir matriz por Entrez (una fila por Entrez con mayor IQR) ---
  # Si el objeto ya trae data_by_entrez, úsalo; si no, constrúyelo desde data_by_gene/anno_by_gene.
  if (!is.null(trans_obj$data_by_entrez) && !is.null(trans_obj$anno_by_entrez)) {
    M_entrez <- as.matrix(trans_obj$data_by_entrez)  # filas=Entrez, cols=muestras
    anno_entrez <- trans_obj$anno_by_entrez
  } else {
    M0   <- as.matrix(trans_obj$data_by_gene)        # filas = best_id, cols = muestras
    anno <- trans_obj$anno_by_gene
    entrez_vec <- as.character(anno$EntrezGeneID)[match(rownames(M0), anno$best_id)]
    keep <- !is.na(entrez_vec) & nzchar(entrez_vec)
    M0   <- M0[keep, , drop = FALSE]
    entrez_vec <- entrez_vec[keep]
    iqr_rows  <- apply(M0, 1, IQR, na.rm = TRUE)
    idx_list  <- split(seq_len(nrow(M0)), entrez_vec)
    pick_rows <- vapply(idx_list, function(ix) ix[which.max(iqr_rows[ix])], integer(1))
    M_entrez  <- M0[pick_rows, , drop = FALSE]
    rownames(M_entrez) <- names(idx_list)            # rownames = Entrez únicos
    anno_entrez <- anno[match(rownames(M_entrez), anno$EntrezGeneID),
                        c("best_id","id_type","GeneSymbol","GeneName",
                          "EnsemblID","EntrezGeneID","RefSeqAccession",
                          "GenbankAccession","Cytoband"), drop = FALSE]
    rownames(anno_entrez) <- anno_entrez$EntrezGeneID
  }
  
  # --- 0b) Mapa de etiquetas Entrez -> GeneSymbol (rellena con Entrez si falta) ---
  id_to_label <- {
    gs <- anno_entrez$GeneSymbol
    names(gs) <- rownames(anno_entrez)              # nombres = Entrez
    gs[is.na(gs) | gs == ""] <- names(gs)[is.na(gs) | gs == ""]
    gs
  }
  
  # --- 0c) Validación/relleno con clusterProfiler::bitr (opcional) ---
  if (isTRUE(validate_with_bitr) && requireNamespace("clusterProfiler", quietly = TRUE)) {
    # Resolver OrgDb: aceptar objeto o nombre de paquete
    OrgDb_obj <- NULL
    if (inherits(orgdb, "OrgDb")) {
      OrgDb_obj <- orgdb
    } else if (is.character(orgdb) && length(orgdb) == 1L) {
      # cargar el paquete si está disponible
      if (!requireNamespace(orgdb, quietly = TRUE)) {
        warning(sprintf("No se pudo cargar '%s'; se omite bitr.", orgdb))
      } else {
        # obtener el objeto OrgDb exportado con el mismo nombre
        OrgDb_obj <- try(getExportedValue(orgdb, orgdb), silent = TRUE)
        if (inherits(OrgDb_obj, "try-error") || !inherits(OrgDb_obj, "OrgDb")) {
          # fallback: buscar cualquier objeto OrgDb exportado
          nms <- getNamespaceExports(orgdb)
          cand <- Filter(function(x) inherits(getExportedValue(orgdb, x), "OrgDb"), nms)
          if (length(cand)) {
            OrgDb_obj <- getExportedValue(orgdb, cand[[1]])
          } else {
            warning(sprintf("No se encontró un objeto OrgDb válido en '%s'; se omite bitr.", orgdb))
            OrgDb_obj <- NULL
          }
        }
      }
    }
    if (!is.null(OrgDb_obj) && inherits(OrgDb_obj, "OrgDb")) {
      # ids a validar (ENTREZ únicamente)
      entrez_ids <- names(id_to_label)
      entrez_ids <- unique(entrez_ids[nzchar(entrez_ids)])
      if (length(entrez_ids)) {
        bit <- try(
          suppressMessages(
            clusterProfiler::bitr(entrez_ids,
                                  fromType = "ENTREZID",
                                  toType   = "SYMBOL",
                                  OrgDb    = OrgDb_obj)
          ),
          silent = TRUE
        )
        if (!inherits(bit, "try-error") && !is.null(bit) && nrow(bit) > 0) {
          map_sym <- setNames(as.character(bit$SYMBOL), as.character(bit$ENTREZID))
          # solo rellenar donde aún el "símbolo" es el propio entrez
          missing <- names(id_to_label)[ id_to_label == names(id_to_label) ]
          repl <- map_sym[missing]
          repl[is.na(repl) | !nzchar(repl)] <- missing
          id_to_label[missing] <- repl
        }
      }
    }
  }
  # --- 1) Pasar a muestras x genes (columnas = Entrez) ---
  Xg <- t(M_entrez)
  mode(Xg) <- "numeric"
  
  # --- 2) Alinear a metadata / y ---
  ids0 <- intersect(rownames(metadata), rownames(Xg))
  if (!length(ids0)) stop("No hay IDs comunes en transcriptómica.")
  metadata <- metadata[ids0, , drop = FALSE]
  Xg <- Xg[ids0, , drop = FALSE]
  if (is.null(names(y))) names(y) <- rownames(metadata)
  y <- droplevels(y[rownames(metadata)])
  if (anyNA(y)) stop("y contiene NA tras alinear.")
  
  # --- 3) Filtrar por IDs fijos (si llegan) ---
  if (!is.null(fixed_keep_ids)) {
    ids_use <- intersect(fixed_keep_ids, rownames(Xg))
    Xg <- Xg[ids_use, , drop = FALSE]
    metadata <- metadata[ids_use, , drop = FALSE]
    y <- droplevels(y[ids_use])
  }
  
  # --- 4) Split + escalado coherente ---
  if (!is.null(fixed_train_ids) && !is.null(fixed_test_ids)) {
    tr <- intersect(fixed_train_ids, rownames(Xg))
    te <- intersect(fixed_test_ids,  rownames(Xg))
    X_tr <- Xg[tr, , drop = FALSE]; X_te <- Xg[te, , drop = FALSE]
    y_tr <- droplevels(y[tr]);       y_te <- droplevels(y[te])
    X_tr_sc <- scale(X_tr)
    ctr <- attr(X_tr_sc, "scaled:center"); scl <- attr(X_tr_sc, "scaled:scale")
    X_te_sc <- scale(X_te, center = ctr, scale = scl)
    ss <- list(mode="split_scaled_fixed",
               X_train=X_tr, X_test=X_te,
               X_train_scaled=X_tr_sc, X_test_scaled=X_te_sc,
               tX_train_scaled=t(X_tr_sc), tX_test_scaled=t(X_te_sc),
               y_train=y_tr, y_test=y_te, center=ctr, scale=scl)
  } else {
    ss <- split_scale(X = Xg, y = y, p = p_train, seed = seed,
                      remove_outliers = FALSE, do_split = TRUE,
                      min_test_per_class = min_test_per_class)
  }
  
  # --- 5) Helper TEST: exige mismas columnas (Entrez) en el mismo orden ---
  transform_test <- function(X_new){
    req_cols <- colnames(ss$X_train_scaled)
    if (!all(req_cols %in% colnames(X_new)))
      stop("X_new no contiene todos los Entrez del train.")
    X_new <- X_new[, req_cols, drop = FALSE]
    X_new <- as.matrix(X_new); mode(X_new) <- "numeric"
    X_new_scaled <- scale(X_new, center = ss$center, scale = ss$scale)
    list(X_new_scaled = X_new_scaled,
         tX_new_scaled = t(X_new_scaled))
  }
  
  list(
    X_all     = Xg,
    split     = ss,
    center    = ss$center,
    scale     = ss$scale,
    genes     = colnames(Xg),     # EntrezGeneID
    metadata  = metadata,
    id_to_label = id_to_label,    # <- mapa Entrez -> símbolo (validado si se pidió)
    transform_test = transform_test
  )
}

process_proteoma <- function(
    path_citoquina,
    path_miokina,
    sep_cito = "",
    sep_mio = "",
    metadata, y,
    fixed_keep_ids = NULL,
    fixed_train_ids = NULL,
    fixed_test_ids  = NULL,
    use_feats_for_outliers = NULL,
    remove_outliers = TRUE,
    p_train = 0.7,
    seed = 123,
    min_test_per_class = 1,
    epsilon = 0
){
  # 1) Carga
  citoquina <- read.delim(path_citoquina, sep = sep_cito, row.names = 1, check.names = FALSE)
  miokina   <- read.delim(path_miokina,   sep = sep_mio,  row.names = 1, check.names = FALSE)
  proteoma  <- cbind(citoquina, miokina)
  
  # 2) Alinear a metadata
  ids0 <- intersect(rownames(metadata), rownames(proteoma))
  if (!length(ids0)) stop("No hay IDs comunes.")
  metadata <- metadata[ids0, , drop = FALSE]
  proteoma <- proteoma[ids0, , drop = FALSE]
  if (is.null(names(y))) names(y) <- rownames(metadata)
  y <- droplevels(y[rownames(metadata)])
  
  # 3) Filtro por IDs fijos o outliers
  if (!is.null(fixed_keep_ids)) {
    ids_use <- intersect(fixed_keep_ids, rownames(proteoma))
    proteoma <- proteoma[ids_use, , drop = FALSE]
    metadata <- metadata[ids_use, , drop = FALSE]
    y        <- droplevels(y[ids_use])
  } else if (isTRUE(remove_outliers)) {
    proteoma <- fun_remove_outliers(proteoma, use_feats_for_outliers)
    ids1 <- intersect(rownames(metadata), rownames(proteoma))
    proteoma <- proteoma[ids1, , drop = FALSE]
    metadata <- metadata[ids1, , drop = FALSE]
    y        <- droplevels(y[ids1])
  }
  if (anyNA(y)) stop("y contiene NA tras filtrar/alinear.")
  
  # 4) Normalización composicional + log2
  rs <- rowSums(proteoma, na.rm = TRUE)
  if (any(rs == 0)) { if (epsilon <= 0) stop("Filas con suma 0."); rs <- rs + epsilon }
  X_log2 <- log2(sweep(proteoma, 1, rs, "/"))
  
  # 5) Split + escalado
  if (!is.null(fixed_train_ids) && !is.null(fixed_test_ids)) {
    tr <- intersect(fixed_train_ids, rownames(X_log2))
    te <- intersect(fixed_test_ids,  rownames(X_log2))
    X_tr <- X_log2[tr, , drop = FALSE]; X_te <- X_log2[te, , drop = FALSE]
    y_tr <- droplevels(y[tr]);           y_te <- droplevels(y[te])
    X_tr_sc <- scale(X_tr)
    ctr <- attr(X_tr_sc, "scaled:center"); scl <- attr(X_tr_sc, "scaled:scale")
    X_te_sc <- scale(X_te, center = ctr, scale = scl)
    ss <- list(mode="split_scaled_fixed",
               X_train=X_tr, X_test=X_te,
               X_train_scaled=X_tr_sc, X_test_scaled=X_te_sc,
               tX_train_scaled=t(X_tr_sc), tX_test_scaled=t(X_te_sc),
               y_train=y_tr, y_test=y_te, center=ctr, scale=scl)
  } else {
    ss <- split_scale(X = X_log2, y = y, p = p_train, seed = seed,
                      remove_outliers = FALSE, do_split = TRUE,
                      min_test_per_class = min_test_per_class)
  }
  
  # 6) Helper para TEST futuro
  transform_test <- function(X_new){
    # mismas columnas y orden que el train
    req_cols <- colnames(ss$X_train_scaled)
    if (!all(req_cols %in% colnames(X_new)))
      stop("X_new no contiene todas las columnas del train.")
    X_new <- X_new[, req_cols, drop = FALSE]
    # misma normalización
    rsn <- rowSums(X_new, na.rm = TRUE)
    if (any(rsn == 0)) { if (epsilon <= 0) stop("X_new filas con suma 0."); rsn <- rsn + epsilon }
    X_new_log2 <- log2(sweep(X_new, 1, rsn, "/"))
    # mismo escalado
    X_new_scaled <- scale(X_new_log2, center = ss$center, scale = ss$scale)
    list(X_new_log2 = X_new_log2,
         X_new_scaled = X_new_scaled,
         tX_new_scaled = t(X_new_scaled))
  }
  
  list(X_log2_all = X_log2,
       split = ss,
       center = ss$center, scale = ss$scale,
       metadata = metadata,
       transform_test = transform_test)
}


split_scale <- function(
    X,                      # matriz/data.frame: n x p (obs x features)
    y = NULL,               # factor/clase para split estratificado (requerido si do_split=TRUE)
    p = 0.7,                # proporción train
    seed = 123,
    remove_outliers = FALSE,
    outlier_patterns = NULL,# vector de patrones regex a buscar en rownames(X)
    outlier_ids = NULL,     # vector exacto de IDs (rownames) a remover
    do_split = NULL,        # si NULL: se infiere (TRUE si remove_outliers, FALSE si no)
    min_test_per_class = 1
) {
  # Coerciones
  rnames <- rownames(X)
  clnames <- colnames(X)
  X <- as.matrix(as.data.frame(lapply(as.data.frame(X),FUN = as.numeric)))
  colnames(X) <- clnames
  rownames(X) <- rnames
  rn <- rownames(X)
  if (is.null(do_split)) do_split <- isTRUE(remove_outliers)
  
  # --- 1) Remoción de outliers (opcional) ---
  to_drop <- logical(nrow(X))
  if (isTRUE(remove_outliers)) {
    if (!is.null(outlier_patterns) && length(outlier_patterns) > 0) {
      if (is.null(rn)) stop("Se requieren rownames(X) para usar outlier_patterns.")
      rx <- paste(outlier_patterns, collapse = "|")
      to_drop <- to_drop | grepl(rx, rn)
    }
    if (!is.null(outlier_ids) && length(outlier_ids) > 0) {
      if (is.null(rn)) stop("Se requieren rownames(X) para usar outlier_ids.")
      to_drop <- to_drop | rn %in% outlier_ids
    }
  }
  removed_ids <- if (!is.null(rn)) rn[to_drop] else as.character(which(to_drop))
  Xf <- X[!to_drop, , drop = FALSE]
  yf <- if (!is.null(y)) y[!to_drop] else NULL
  
  # --- 2) Sin split: solo escalar todo con sus propios parámetros ---
  if (!isTRUE(do_split)) {
    X_scaled <- scale(Xf)
    ctr <- attr(X_scaled, "scaled:center")
    scl <- attr(X_scaled, "scaled:scale")
    tX_scaled <- t(X_scaled) # "traspuesta de X" ya escalada por columnas originales
    
    return(list(
      mode = "nosplit_scaled",
      removed_ids = removed_ids,
      X_scaled = X_scaled,
      tX_scaled = tX_scaled,
      center = ctr,
      scale = scl,
      n = nrow(Xf),
      p = ncol(Xf)
    ))
  }
  
  # --- 3) Con split estratificado ---
  if (is.null(yf)) stop("y es requerido para split estratificado.")
  if (anyNA(yf)) stop("y contiene NA.")
  yf <- as.factor(as.character(yf))

  set.seed(seed)
  idx_train <- integer(0)
  for (lev in levels(yf)) {
    idx <- which(yf == lev)
    n <- length(idx)
    if (n <= min_test_per_class) {
      stop(sprintf("Clase '%s' tiene %d obs (< min_test_per_class)", lev, n))
    }
    n_train <- max(1, min(n - min_test_per_class, round(p * n)))
    idx_train <- c(idx_train, sample(idx, n_train, replace = FALSE))
  }
  idx_train <- sort(idx_train)
  idx_test  <- setdiff(seq_len(nrow(Xf)), idx_train)
  
  X_train <- Xf[idx_train, , drop = FALSE]
  X_test  <- Xf[idx_test,  , drop = FALSE]
  y_train <- yf[idx_train]
  y_test  <- yf[idx_test]
  
  # --- 4) Escalado: train con sus propios parámetros, test con atributos de train ---
  X_train_sc <- scale(X_train)
  ctr <- attr(X_train_sc, "scaled:center")
  scl <- attr(X_train_sc, "scaled:scale")
  
  X_test_sc <- scale(X_test, center = ctr, scale = scl)
  
  # Traspuestas escaladas coherentes: usar las versiones ya escaladas y luego transponer
  tX_train_sc <- t(X_train_sc)
  tX_test_sc  <- t(X_test_sc)
  
  # --- 5) Salida ---
  list(
    mode = "split_scaled",
    removed_ids = removed_ids,
    idx_train = idx_train,
    idx_test = idx_test,
    X_train = X_train,           # sin escalar
    X_test  = X_test,            # sin escalar
    X_train_scaled = X_train_sc, # escalado con center/scale de train
    X_test_scaled  = X_test_sc,  # escalado con center/scale de train
    tX_train_scaled = tX_train_sc,
    tX_test_scaled  = tX_test_sc,
    y_train = y_train,
    y_test  = y_test,
    counts = list(
      total = table(yf),
      train = table(y_train),
      test  = table(y_test)
    ),
    center = ctr,
    scale  = scl,
    n_train = nrow(X_train),
    n_test  = nrow(X_test),
    p = ncol(Xf)
  )
}




# Congela IDs y split global para todos los bloques
define_split_ids <- function(
    samples_metadata,            # data.frame con rownames = IDs y columna de grupo
    grupo_col = "grupo",
    outlier_patterns = c("NP_02","NP_03","SO_14","6M"),
    outlier_ids = NULL,
    p = 0.6,
    seed = 42,
    min_test_per_class = 1
){
  if (is.null(rownames(samples_metadata)))
    stop("samples_metadata necesita rownames=IDs de muestra.")
  if (!grupo_col %in% colnames(samples_metadata))
    stop(sprintf("No existe la columna '%s' en samples_metadata.", grupo_col))
  
  ids_master <- rownames(samples_metadata)
  
  # y nombrado
  y_master <- setNames(as.factor(as.character(samples_metadata[[grupo_col]])), ids_master)
  if (anyNA(y_master)) stop("y_master contiene NA.")
  
  # outliers por patrón + lista explícita
  rm_pat <- if (length(outlier_patterns) > 0)
    ids_master[grepl(paste(outlier_patterns, collapse="|"), ids_master)] else character(0)
  rm_ids <- unique(c(rm_pat, intersect(ids_master, outlier_ids %||% character(0))))
  
  keep_ids <- setdiff(ids_master, rm_ids)
  if (length(keep_ids) == 0) stop("No quedan IDs tras remover outliers.")
  
  # split estratificado
  set.seed(seed)
  y_keep <- droplevels(y_master[keep_ids])
  
  idx_tr <- integer(0)
  for (lev in levels(y_keep)) {
    idx <- which(y_keep == lev)
    n <- length(idx)
    if (n <= min_test_per_class)
      stop(sprintf("Clase '%s' tiene %d obs (< min_test_per_class).", lev, n))
    n_tr <- max(1, min(n - min_test_per_class, round(p * n)))
    idx_tr <- c(idx_tr, sample(idx, n_tr, replace = FALSE))
  }
  idx_tr <- sort(idx_tr)
  idx_te <- setdiff(seq_along(keep_ids), idx_tr)
  
  train_ids <- keep_ids[idx_tr]
  test_ids  <- keep_ids[idx_te]
  
  list(
    ids_master = ids_master,
    y_master   = y_master,
    removed_ids = rm_ids,
    keep_ids   = keep_ids,
    train_ids  = train_ids,
    test_ids   = test_ids,
    p = p, seed = seed
  )
}

# helper para operador coalescencia
`%||%` <- function(a,b) if (!is.null(a)) a else b

process_metaboloma <- function(
    path_pos, path_neg,
    metadata, y,
    sep_pos = "", sep_neg = "",
    prefix_pos = "pos_", prefix_neg = "neg_",
    fixed_keep_ids = NULL,
    fixed_train_ids = NULL,
    fixed_test_ids  = NULL,
    align_by_suffix = FALSE,                 # TRUE si necesitas casar por sufijo numérico
    suffix_pattern  = ".*_(\\d+)",           # como en tu script
    p_train = 0.7, seed = 123, min_test_per_class = 1
){
  # --- leer POS ---
  df_raw <- readxl::read_excel(path_pos, col_names = FALSE)
  sample_names <- as.character(df_raw[1, -1])
  data_matrix <- as.data.frame(df_raw[-c(1, 2), ])
  rownames(data_matrix) <- data_matrix[[1]]
  data_matrix <- data_matrix[, -1, drop = FALSE]
  colnames(data_matrix) <- sample_names
  data_matrix[] <- lapply(data_matrix, as.numeric)
  df_final <- as.data.frame(t(data_matrix))
  colnames(df_final) <- paste0(prefix_pos, colnames(df_final))
  
  # --- leer NEG ---
  df_raw2 <- readxl::read_excel(path_neg, col_names = FALSE)
  sample_names2 <- as.character(df_raw2[1, -1])
  data_matrix2 <- as.data.frame(df_raw2[-c(1, 2), ])
  rownames(data_matrix2) <- data_matrix2[[1]]
  data_matrix2 <- data_matrix2[, -1, drop = FALSE]
  colnames(data_matrix2) <- sample_names2
  data_matrix2[] <- lapply(data_matrix2, as.numeric)
  df_final2 <- as.data.frame(t(data_matrix2))
  colnames(df_final2) <- paste0(prefix_neg, colnames(df_final2))
  
  # --- combinar ---
  metaboloma <- cbind(df_final, df_final2)
  if (anyNA(metaboloma)) {
    keep_cols <- which(colSums(is.na(metaboloma)) < nrow(metaboloma))
    metaboloma <- metaboloma[, keep_cols, drop = FALSE]
  }
  rownames(metaboloma) <- gsub("-", "_", rownames(df_final))
  
  # --- alinear a metadata / IDs fijos ---
  if (is.null(names(y))) names(y) <- rownames(metadata)
  if (!align_by_suffix) {
    ids0 <- intersect(rownames(metadata), rownames(metaboloma))
    if (!length(ids0)) stop("No hay IDs comunes entre metadata y metaboloma.")
    metadata <- metadata[ids0, , drop = FALSE]
    metaboloma <- metaboloma[ids0, , drop = FALSE]
    y <- droplevels(y[ids0])
  } else {
    ref_ids <- if (!is.null(fixed_keep_ids)) fixed_keep_ids else rownames(metadata)
    suf_ref <- as.numeric(gsub(suffix_pattern, "\\1", ref_ids))
    suf_met <- as.numeric(gsub(suffix_pattern, "\\1", rownames(metaboloma)))
    idx <- match(suf_ref, suf_met)
    if (anyNA(idx)) stop("No se pudieron casar sufijos entre metaboloma y referencia.")
    metaboloma <- metaboloma[idx, , drop = FALSE]
    rownames(metaboloma) <- ref_ids
    metadata <- metadata[ref_ids, , drop = FALSE]
    y <- droplevels(y[ref_ids])
  }
  if (anyNA(y)) stop("y contiene NA tras la alineación.")
  
  # --- filtrar por IDs fijos si se dan ---
  if (!is.null(fixed_keep_ids)) {
    ids_use <- intersect(fixed_keep_ids, rownames(metaboloma))
    metaboloma <- metaboloma[ids_use, , drop = FALSE]
    metadata   <- metadata[ids_use, , drop = FALSE]
    y          <- droplevels(y[ids_use])
  }
  
  # --- split + escalado coherente con otros bloques ---
  if (!is.null(fixed_train_ids) && !is.null(fixed_test_ids)) {
    tr <- intersect(fixed_train_ids, rownames(metaboloma))
    te <- intersect(fixed_test_ids,  rownames(metaboloma))
    X_tr <- metaboloma[tr, , drop = FALSE]; X_te <- metaboloma[te, , drop = FALSE]
    y_tr <- droplevels(y[tr]);              y_te <- droplevels(y[te])
    X_tr_sc <- scale(X_tr)
    ctr <- attr(X_tr_sc, "scaled:center"); scl <- attr(X_tr_sc, "scaled:scale")
    X_te_sc <- scale(X_te, center = ctr, scale = scl)
    ss <- list(
      mode="split_scaled_fixed",
      X_train=X_tr, X_test=X_te,
      X_train_scaled=X_tr_sc, X_test_scaled=X_te_sc,
      tX_train_scaled=t(X_tr_sc), tX_test_scaled=t(X_te_sc),
      y_train=y_tr, y_test=y_te,
      center=ctr, scale=scl
    )
  } else {
    ss <- split_scale(
      X = metaboloma, y = y, p = p_train, seed = seed,
      remove_outliers = FALSE, do_split = TRUE,
      min_test_per_class = min_test_per_class
    )
  }
  
  # --- helper para TEST futuro ---
  transform_test <- function(X_new){
    if (!all(colnames(ss$X_train) %in% colnames(X_new)))
      stop("X_new no contiene todas las variables del train metabolómico.")
    X_new <- X_new[, colnames(ss$X_train), drop = FALSE]
    X_new <- as.matrix(X_new); mode(X_new) <- "numeric"
    X_new_scaled <- scale(X_new, center = ss$center, scale = ss$scale)
    list(
      X_new_scaled   = X_new_scaled,
      tX_new_scaled  = t(X_new_scaled)
    )
  }
  
  list(
    X_all = metaboloma,
    split = ss,
    center = ss$center,
    scale  = ss$scale,
    metadata = metadata,
    transform_test = transform_test
  )
}

process_clinicos <- function(
    path_excel,
    metadata, y,
    fixed_keep_ids = NULL,
    fixed_train_ids = NULL,
    fixed_test_ids  = NULL,
    p_train = 0.7, seed = 123, min_test_per_class = 1,
    vars_clinicas = NULL
){
  raw <- readxl::read_excel(path_excel, col_names = FALSE)
  
  # -----------------------------
  # 1) Subset y nombres
  # -----------------------------
  dat <- raw[4:32, 4:44]
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  
  nombres_sujetos <- as.character(raw[4:32, 3][[1]])
  col_names <- as.character(unlist(raw[3, 4:44]))
  
  rownames(dat) <- nombres_sujetos
  names(dat) <- col_names
  
  # -----------------------------
  # 2) Conversión de tipos
  # -----------------------------
  num_cols <- setdiff(names(dat), c("Grupo", "Género"))
  dat[num_cols] <- lapply(dat[num_cols], as.numeric)
  
  dat <- dplyr::mutate(
    dat,
    Grupo  = factor(Grupo),
    Género = factor(Género)
  )
  # -----------------------------
  # 3) Corrección original filas 20:22
  # -----------------------------
  idx_fix <- intersect(20:22, seq_len(nrow(dat)))
  
  if (length(idx_fix) > 0) {
    cols_fix <- intersect(c("Grupo", "Género"), colnames(dat))
    
    dat[idx_fix, cols_fix] <-
      lapply(dat[idx_fix, cols_fix, drop = FALSE], function(x) {
        factor(replace(as.character(x), x == "sobrepeso", "obeso"))
      })
  }
  
  # -----------------------------
  # 4) Variables clínicas seleccionadas
  # -----------------------------
  if (is.null(vars_clinicas)) {
    stop("Debes pasar vars_clinicas desde la configuración del análisis.")
  }
  
  vars <- vars_clinicas
  vars_missing <- setdiff(vars, colnames(dat))
  
  if (length(vars_missing) > 0) {
    stop(
      "Faltan variables clínicas en el Excel: ",
      paste(vars_missing, collapse = ", ")
    )
  }
  
  dat_final <- as.data.frame(dat[, vars, drop = FALSE])
  rownames(dat_final) <- nombres_sujetos
  
  # -----------------------------
  # 5) Alineación con metadata/y
  # -----------------------------
  ids0 <- intersect(rownames(metadata), rownames(dat_final))
  
  if (!length(ids0)) {
    stop("No hay IDs comunes entre metadata y clínicos.")
  }
  
  metadata  <- metadata[ids0, , drop = FALSE]
  dat_final <- dat_final[ids0, , drop = FALSE]
  
  if (is.null(names(y))) {
    names(y) <- rownames(metadata)
  }
  
  y <- droplevels(y[rownames(metadata)])
  
  if (anyNA(y)) {
    stop("y contiene NA tras alinear clínicos.")
  }
  
  # -----------------------------
  # 6) Matriz numérica clínica
  # -----------------------------
  num_df <- as.data.frame(
    lapply(dat_final, as.numeric),
    check.names = FALSE
  )
  
  rownames(num_df) <- rownames(dat_final)
  
  # No imputar: si hay NA/Inf, parar
  X_raw <- as.matrix(num_df)
  mode(X_raw) <- "numeric"
  
  if (anyNA(X_raw)) {
    
    na_pos <- which(is.na(X_raw), arr.ind = TRUE)
    
    na_tbl <- data.frame(
      sample = rownames(X_raw)[na_pos[, "row"]],
      variable = colnames(X_raw)[na_pos[, "col"]],
      value_original = as.character(dat_final[cbind(na_pos[, "row"], na_pos[, "col"])]),
      stringsAsFactors = FALSE
    )
    
    cat("\n============================================================\n")
    cat("ERROR: NA en clínicos numéricos después de as.numeric()\n")
    cat("============================================================\n")
    print(na_tbl, row.names = FALSE)
    
    cat("\nResumen por variable:\n")
    print(sort(table(na_tbl$variable), decreasing = TRUE))
    
    cat("\nResumen por muestra:\n")
    print(sort(table(na_tbl$sample), decreasing = TRUE))
    
    write.csv(
      na_tbl,
      file.path(OUTPUT_DIR, "clinicos_NA_detectados.csv"),
      row.names = FALSE
    )
    
    stop(
      "Hay NA en clínicos numéricos. ",
      "Revisa el archivo: ",
      file.path(OUTPUT_DIR, "clinicos_NA_detectados.csv")
    )
  }
  
  if (any(!is.finite(X_raw))) {
    stop("Hay Inf/NaN en clínicos numéricos. No se imputa nada.")
  }
  
  # -----------------------------
  # 7) Keep IDs fijos
  # -----------------------------
  if (!is.null(fixed_keep_ids)) {
    ids_use <- intersect(fixed_keep_ids, rownames(X_raw))
    
    if (!length(ids_use)) {
      stop("fixed_keep_ids no coincide con ninguna muestra clínica.")
    }
    
    X_raw    <- X_raw[ids_use, , drop = FALSE]
    metadata <- metadata[ids_use, , drop = FALSE]
    y        <- droplevels(y[ids_use])
  }
  
  if (anyNA(y)) {
    stop("y contiene NA tras fixed_keep_ids.")
  }
  
  # -----------------------------
  # 8) Transformación seleccionada:
  #    Normal + log natural
  #    Escalado SIN fuga:
  #    - center/scale calculados SOLO en TRAIN
  #    - TEST escalado con parámetros de TRAIN
  # -----------------------------
  if (any(X_raw <= 0)) {
    vars_bad <- colnames(X_raw)[colSums(X_raw <= 0) > 0]
    
    stop(
      "Hay valores <= 0 en clínicos numéricos. ",
      "No se puede aplicar log() directo en: ",
      paste(vars_bad, collapse = ", ")
    )
  }
  
  X_log <- log(X_raw)
  
  if (anyNA(X_log) || any(!is.finite(X_log))) {
    stop("X_log contiene NA/Inf/NaN tras log().")
  }
  
  # -----------------------------
  # 9) Split: respeta IDs fijos
  # -----------------------------
  if (!is.null(fixed_train_ids) && !is.null(fixed_test_ids)) {
    
    tr <- fixed_train_ids[fixed_train_ids %in% rownames(X_log)]
    te <- fixed_test_ids [fixed_test_ids  %in% rownames(X_log)]
    
    if (!length(tr) || !length(te)) {
      stop("IDs fijos no encontrados en clínicos.")
    }
    
    X_tr <- X_log[tr, , drop = FALSE]
    X_te <- X_log[te, , drop = FALSE]
    
    y_tr <- droplevels(y[tr])
    y_te <- droplevels(y[te])
    
  } else {
    
    set.seed(seed)
    y_fac <- droplevels(y)
    
    idx_train <- integer(0)
    
    for (lev in levels(y_fac)) {
      idx <- which(y_fac == lev)
      n <- length(idx)
      
      if (n <= min_test_per_class) {
        stop(sprintf("Clase '%s' tiene %d obs (< min_test_per_class)", lev, n))
      }
      
      n_tr <- max(1, min(n - min_test_per_class, round(p_train * n)))
      idx_train <- c(idx_train, sample(idx, n_tr, replace = FALSE))
    }
    
    idx_train <- sort(idx_train)
    idx_test  <- setdiff(seq_len(nrow(X_log)), idx_train)
    
    X_tr <- X_log[idx_train, , drop = FALSE]
    X_te <- X_log[idx_test,  , drop = FALSE]
    
    y_tr <- y_fac[idx_train]
    y_te <- y_fac[idx_test]
    
    tr <- rownames(X_tr)
    te <- rownames(X_te)
  }
  
  # -----------------------------
  # 10) Center/scale SOLO con TRAIN
  # -----------------------------
  sds_tr <- apply(X_tr, 2, sd)
  
  if (any(sds_tr == 0 | !is.finite(sds_tr))) {
    vars_const <- names(sds_tr)[sds_tr == 0 | !is.finite(sds_tr)]
    
    stop(
      "Hay variables clínicas con varianza 0 o SD no finita en TRAIN tras log(): ",
      paste(vars_const, collapse = ", ")
    )
  }
  
  X_tr_sc <- scale(X_tr, center = TRUE, scale = TRUE)
  
  ctr <- attr(X_tr_sc, "scaled:center")
  scl <- attr(X_tr_sc, "scaled:scale")
  
  X_te_sc <- scale(X_te, center = ctr, scale = scl)
  
  # Matriz ALL escalada con parámetros de TRAIN.
  # Esto es útil para diagnóstico, pero center/scale vienen solo del TRAIN.
  X_scaled_all <- scale(X_log, center = ctr, scale = scl)
  
  if (anyNA(X_tr_sc) || any(!is.finite(X_tr_sc))) {
    stop("X_train_scaled contiene NA/Inf/NaN.")
  }
  
  if (anyNA(X_te_sc) || any(!is.finite(X_te_sc))) {
    stop("X_test_scaled contiene NA/Inf/NaN.")
  }
  
  if (anyNA(X_scaled_all) || any(!is.finite(X_scaled_all))) {
    stop("X_scaled_all contiene NA/Inf/NaN.")
  }
  
  ss <- list(
    mode = if (!is.null(fixed_train_ids) && !is.null(fixed_test_ids)) {
      "split_log_scaled_fixed_train_only"
    } else {
      "split_log_scaled_train_only"
    },
    X_train = X_tr,
    X_test  = X_te,
    X_train_scaled = X_tr_sc,
    X_test_scaled  = X_te_sc,
    tX_train_scaled = t(X_tr_sc),
    tX_test_scaled  = t(X_te_sc),
    y_train = y_tr,
    y_test  = y_te,
    center = ctr,
    scale  = scl
  )
  
  if (exists("idx_train")) {
    ss$idx_train <- idx_train
    ss$idx_test  <- idx_test
  }
  # -----------------------------
  # 10) Transformación para datos nuevos
  # -----------------------------
  transform_test <- function(X_new) {
    
    req_cols <- colnames(X_log)
    
    if (!all(req_cols %in% colnames(X_new))) {
      faltan <- setdiff(req_cols, colnames(X_new))
      stop(
        "X_new no contiene todas las columnas clínicas requeridas: ",
        paste(faltan, collapse = ", ")
      )
    }
    
    X_new <- X_new[, req_cols, drop = FALSE]
    X_new <- as.matrix(X_new)
    mode(X_new) <- "numeric"
    
    if (anyNA(X_new)) {
      stop("X_new contiene NA. No se imputa nada.")
    }
    
    if (any(!is.finite(X_new))) {
      stop("X_new contiene Inf/NaN. No se imputa nada.")
    }
    
    if (any(X_new <= 0)) {
      vars_bad <- colnames(X_new)[colSums(X_new <= 0) > 0]
      stop(
        "X_new contiene valores <= 0. No se puede aplicar log() en: ",
        paste(vars_bad, collapse = ", ")
      )
    }
    
    X_new_log <- log(X_new)
    X_new_scaled <- scale(X_new_log, center = ctr, scale = scl)
    
    list(
      X_new_log = X_new_log,
      X_new_scaled = X_new_scaled,
      tX_new_scaled = t(X_new_scaled)
    )
  }  
  # -----------------------------
  # 11) Salida
  # -----------------------------
  list(
    X_raw_all    = X_raw,
    X_log_all    = X_log,
    X_scaled_all = X_scaled_all,
    
    # Alias para compatibilidad con funciones antiguas
    # OJO: ahora X_log_all sí es log sin escalar
    split     = ss,
    center    = ctr,
    scale     = scl,
    metadata  = metadata,
    cols_numeric = colnames(X_log),
    transform_used = "log_natural_plus_center_scale",
    transform_test = transform_test
  )
}

save_block_plots <- function(
    res, out_dir, name,
    group_col = "grupo",
    label_col = "CÓDIGO.FINAL",
    palette = c(FD="red", SO="blue", NP="black"),
    width = 7, height = 5, dpi = 300, quality = 95,
    print_plots = TRUE
){
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  if ("X_log2_all" %in% names(res)) X0 <- res$X_log2_all
  else if ("X_log_all" %in% names(res)) X0 <- res$X_log_all
  else if ("X_all" %in% names(res)) X0 <- res$X_all
  else stop("El objeto no contiene X_log2_all, X_log_all ni X_all.")
  
  X0 <- as.matrix(X0)
  meta <- res$metadata
  stopifnot(all(rownames(X0) %in% rownames(meta)))
  meta_all <- meta[rownames(X0), , drop = FALSE]
  if (!group_col %in% colnames(meta_all)) stop("group_col no existe en metadata.")
  
  ctr <- res$center; scl <- res$scale
  if (is.null(ctr) || is.null(scl)) stop("Faltan center/scale en res.")
  Xz_all <- scale(X0[, names(ctr), drop = FALSE], center = ctr, scale = scl)
  
  tr_ids <- rownames(res$split$X_train)
  X0_tr  <- X0[tr_ids, , drop = FALSE]
  Xz_tr  <- scale(X0_tr[, names(ctr), drop = FALSE], center = ctr, scale = scl)
  meta_tr <- meta[tr_ids, , drop = FALSE]
  
  # => filename primero, luego subset/proc
  .save <- function(p, fn, subset_type, proc_type){
    dir_target <- file.path(out_dir, name, subset_type, proc_type)
    if (!dir.exists(dir_target)) dir.create(dir_target, recursive = TRUE)
    ggsave(
      filename = fn, path = dir_target, plot = p,
      width = width, height = height, units = "in",
      device = "jpeg", dpi = dpi, quality = quality
    )
    if (isTRUE(print_plots)) print(p)
  }
  
  .dens <- function(X, meta_df, ttl){
    df <- cbind(meta_df[, group_col, drop = FALSE], as.data.frame(X))
    colnames(df)[1] <- "grupo"
    mlt <- reshape2::melt(df, id.vars = "grupo")
    ggplot(mlt, aes(value, color = grupo)) +
      geom_density(aes(y = after_stat(density)), alpha = 0.25, linewidth = 0.7) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")) +
      labs(title = ttl, x = "Valor", y = "Densidad", color = "Grupo") +
      scale_color_manual(values = palette)
  }
  
  .pca <- function(X, meta_df, ttl){
    pc <- prcomp(X, center = FALSE, scale. = FALSE)
    var <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 2)
    lbl <- if (label_col %in% colnames(meta_df)) meta_df[[label_col]] else rownames(meta_df)
    plt <- cbind(as.data.frame(pc$x[, 1:2, drop=FALSE]),
                 grupo = as.factor(meta_df[[group_col]]),
                 label = lbl)
    ggplot(plt, aes(PC1, PC2, color = grupo, label = label)) +
      geom_point(size = 2, alpha = 0.9) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")) +
      labs(title = ttl,
           x = paste0("PC1: ", var[1], "%"),
           y = paste0("PC2: ", var[2], "%"),
           color = "Grupo") +
      scale_color_manual(values = palette)
  }
  
  .plsda <- function(X, meta_df, ttl){
    if (!requireNamespace("mixOmics", quietly = TRUE)) stop("mixOmics no instalado.")
    y <- as.factor(meta_df[[group_col]])
    mdl <- mixOmics::plsda(X, y, ncomp = 2, scale = FALSE)
    var <- round(100 * mdl$prop_expl_var$X, 2)
    sc <- as.data.frame(mdl$variates$X)
    sc$grupo <- y
    sc$label <- if (label_col %in% colnames(meta_df)) meta_df[[label_col]] else rownames(meta_df)
    ggplot(sc, aes(comp1, comp2, color = grupo, label = label)) +
      geom_point(size = 2, alpha = 0.9) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")) +
      labs(title = ttl,
           x = paste0("PC1: ", var[1], "%"),
           y = paste0("PC2: ", var[2], "%"),
           color = "Grupo") +
      scale_color_manual(values = palette)
  }
  
  # ALL
  .save(.dens(X0,     meta_all, paste0(name, " · Densidad · ALL · sin procesar")),  "density_all.jpg",  "ALL",  "raw")
  .save(.dens(Xz_all, meta_all, paste0(name, " · Densidad · ALL · procesado")),     "density_all.jpg",  "ALL",  "processed")
  .save(.pca(X0,      meta_all, paste0(name, " · PCA · ALL · sin procesar")),       "pca_all.jpg",      "ALL",  "raw")
  .save(.pca(Xz_all,  meta_all, paste0(name, " · PCA · ALL · procesado")),          "pca_all.jpg",      "ALL",  "processed")
  .save(.plsda(X0,    meta_all, paste0(name, " · PLS-DA · ALL · sin procesar")),    "plsda_all.jpg",    "ALL",  "raw")
  .save(.plsda(Xz_all,meta_all, paste0(name, " · PLS-DA · ALL · procesado")),       "plsda_all.jpg",    "ALL",  "processed")
  
  # TRAIN
  .save(.dens(X0_tr,   meta_tr, paste0(name, " · Densidad · TRAIN · sin procesar")), "density_train.jpg","TRAIN","raw")
  .save(.dens(Xz_tr,   meta_tr, paste0(name, " · Densidad · TRAIN · procesado")),    "density_train.jpg","TRAIN","processed")
  .save(.pca(X0_tr,    meta_tr, paste0(name, " · PCA · TRAIN · sin procesar")),      "pca_train.jpg",    "TRAIN","raw")
  .save(.pca(Xz_tr,    meta_tr, paste0(name, " · PCA · TRAIN · procesado")),         "pca_train.jpg",    "TRAIN","processed")
  .save(.plsda(X0_tr,  meta_tr, paste0(name, " · PLS-DA · TRAIN · sin procesar")),   "plsda_train.jpg",  "TRAIN","raw")
  .save(.plsda(Xz_tr,  meta_tr, paste0(name, " · PLS-DA · TRAIN · procesado")),      "plsda_train.jpg",  "TRAIN","processed")
}


# --- 1) PLS-DA SOLO CON TRAIN: loadings + histogramas (comp1 y comp2) ---
run_plsda_importance_train <- function(
    res, name, out_dir = "plots",
    id_to_label = NULL,
    palette = c(FD="red", SO="blue", NP="black")
){
  if (!requireNamespace("mixOmics", quietly = TRUE))
    stop("mixOmics no instalado.")
  Xtr <- res$split$X_train_scaled
  ytr <- res$split$y_train
  if (is.null(Xtr) || is.null(ytr)) stop("Falta X_train_scaled o y_train.")
  
  mdl <- mixOmics::plsda(Xtr, ytr, ncomp = 2, scale = FALSE)
  
  L_raw <- data.frame(
    variable = rownames(mdl$loadings$X),
    comp1 = as.numeric(mdl$loadings$X[,1]),
    comp2 = as.numeric(mdl$loadings$X[,2]),
    stringsAsFactors = FALSE
  )
  
  get_group_from_plotLoadings <- function(m, comp) {
    aux <- mixOmics::plotLoadings(m, comp = comp, contrib = "max", plot = FALSE)
    mat <- NULL
    if (!is.null(aux$mat) && is.matrix(aux$mat)) mat <- aux$mat
    else if (!is.null(aux$X) && is.data.frame(aux$X)) {
      bool_cols <- vapply(aux$X, function(z) is.logical(z) || all(z %in% c(0,1,TRUE,FALSE), na.rm=TRUE), TRUE)
      if (any(bool_cols)) mat <- as.matrix(aux$X[, bool_cols, drop=FALSE])
    }
    if (is.null(mat)) stop("plotLoadings no devolvió matriz lógica de clases.")
    if (is.null(rownames(mat))) rownames(mat) <- rownames(m$loadings$X)
    grp <- colnames(mat)[max.col(mat, ties.method = "first")]
    data.frame(variable = rownames(mat), comp = comp, grupo_raw = grp, row.names = NULL)
  }
  G <- rbind(get_group_from_plotLoadings(mdl, 1), get_group_from_plotLoadings(mdl, 2))
  
  # limpiar nombres de grupo: "Contrib.FD" -> "FD"
  clean_group <- function(x){
    x <- as.character(x)
    x <- gsub("^(?i)(contrib\\.|class\\.|group\\.)", "", x, perl = TRUE)
    toupper(trimws(x))
  }
  G$grupo <- clean_group(G$grupo_raw)
  
  win <- L_raw |>
    dplyr::mutate(abs1 = abs(comp1), abs2 = abs(comp2)) |>
    dplyr::transmute(variable,
                     comp_win = ifelse(abs1 >= abs2, 1L, 2L),
                     loading_win = ifelse(abs1 >= abs2, comp1, comp2),
                     abs_loading_win = pmax(abs1, abs2)) |>
    dplyr::left_join(G[, c("variable","comp","grupo")], by = c("variable","comp_win" = "comp"))
  
  # Etiquetas legibles
  if (!is.null(id_to_label)) {
    lab <- id_to_label[win$variable]
    lab[is.na(lab) | lab == ""] <- win$variable
    win$label_var <- lab
  } else {
    win$label_var <- win$variable
  }
  
  # Importancia (para la siguiente fase)
  importance_train <- setNames(win$abs_loading_win, win$variable)
  
  # Colores por grupo (FD/SO/NP)
  grupos_presentes <- unique(win$grupo)
  pal_use <- palette[names(palette) %in% grupos_presentes]
  # si aparece algún grupo “raro”, color gris
  extra <- setdiff(grupos_presentes, names(palette))
  if (length(extra) > 0) pal_use <- c(pal_use, setNames(rep("grey70", length(extra)), extra))
  
  # Top-40
  top40 <- win |>
    dplyr::arrange(dplyr::desc(abs_loading_win)) |>
    dplyr::slice(1:40) |>
    dplyr::arrange(loading_win)
  
  p_top40 <- ggplot(top40,
                    aes(x = loading_win,
                        y = factor(label_var, levels = top40$label_var),
                        fill = factor(grupo, levels = names(pal_use)))) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = pal_use, breaks = names(pal_use), labels = names(pal_use), drop = FALSE) +
    theme_minimal(base_size = 11) +
    labs(title = paste0("Top-40 loadings (con signo) · TRAIN · ", name),
         x = "Loading", y = NULL, fill = "Grupo") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())
  
  dir_base <- file.path(out_dir, name, "feature_prefilter_train")
  if (!dir.exists(dir_base)) dir.create(dir_base, recursive = TRUE)
  f_top40 <- file.path(dir_base, "top40_importancia_piramidal_signo.png")
  ggsave(f_top40, p_top40, width = 8, height = 8, dpi = 300)
  print(p_top40)
  
  invisible(list(
    model_train = mdl,
    importance_train = importance_train,
    top40_table = top40[, c("variable","label_var","loading_win","abs_loading_win","grupo","comp_win")],
    paths = list(top40 = f_top40)
  ))
}



# --- helper para obtener la matriz BASE "normalizada sin escalar" por ómica ---
.get_base_matrix_unscaled <- function(res) {
  
  if ("X_log2_all" %in% names(res)) {
    
    return(as.matrix(res$X_log2_all))
    
  } else if ("X_log_all" %in% names(res)) {
    
    return(as.matrix(res$X_log_all))
    
  } else if ("X_all" %in% names(res)) {
    
    return(as.matrix(res$X_all))
    
  } else {
    
    stop("No encuentro X_log2_all / X_log_all / X_all en 'res'.")
  }
}
# --- 2) UMBRAL SOLO CON TRAIN, re-normaliza por ómica, aplica a TEST, refit y plots ---
# === Extensión: devolver también muestras×vars y un "paquete final" coherente ===
apply_importance_threshold_and_refit_from_train <- function(
    res, importance_train, q = 0.9, name, out_dir = "plots") {
  
  if (!requireNamespace("mixOmics", quietly = TRUE))
    stop("mixOmics no instalado.")
  if (is.null(names(importance_train)))
    stop("importance_train debe estar nombrado por variable (nombres = variables).")
  
  # Umbral y selección SOLO usando TRAIN
  thr <- as.numeric(stats::quantile(importance_train, q))
  sel_vars <- names(importance_train)[importance_train >= thr]
  if (!length(sel_vars)) stop("Umbral demasiado alto: 0 variables seleccionadas.")
  
  # Matriz base sin escalar (según ómica)
  X_base <- .get_base_matrix_unscaled(res)            # obs x vars
  keep <- intersect(colnames(X_base), sel_vars)
  if (!length(keep)) stop("Las variables seleccionadas no están en la matriz base.")
  X_base_sel <- X_base[, keep, drop = FALSE]
  
  # IDs
  tr_ids <- rownames(res$split$X_train_scaled)
  te_ids <- rownames(res$split$X_test_scaled)
  if (is.null(tr_ids) || is.null(te_ids))
    stop("Faltan IDs de train/test en res.")
  
  # Re-escalado NUEVO con SOLO TRAIN
  X_tr_base <- X_base_sel[tr_ids, , drop = FALSE]
  X_tr_sc   <- scale(X_tr_base)
  ctr_new   <- attr(X_tr_sc, "scaled:center")
  scl_new   <- attr(X_tr_sc, "scaled:scale")
  X_all_sc  <- scale(X_base_sel, center = ctr_new, scale = scl_new)
  X_te_sc   <- X_all_sc[te_ids, , drop = FALSE]
  
  # Refit PLS-DA SOLO en TRAIN filtrado
  y_tr <- res$split$y_train
  mdl2 <- mixOmics::plsda(X_tr_sc, y_tr, ncomp = 2, scale = FALSE)
  var_pc <- round(100 * mdl2$prop_expl_var$X, 2)
  
  # ---- Salidas en ambas orientaciones ----
  # 1) vars×muestras (como ya usabas para gráficos/diagnóstico)
  X_vs_train <- t(X_tr_sc)
  X_vs_test  <- t(X_te_sc)
  
  # 2) muestras×vars (para detectar par y alimentar modelos)
  X_sv_train <- X_tr_sc
  X_sv_test  <- X_te_sc
  stopifnot(identical(colnames(X_sv_train), colnames(X_sv_test)))
  stopifnot(length(intersect(rownames(X_sv_train), rownames(X_sv_test))) == 0)
  
  # ---- Plots de umbral y scores (opcional, igual que antes) ----
  dir_base <- file.path(out_dir, name, "feature_filtering_from_train")
  if (!dir.exists(dir_base)) dir.create(dir_base, recursive = TRUE)
  
  df_imp <- data.frame(variable = names(importance_train),
                       abs_loading = as.numeric(importance_train))
  p_hist_thr <- ggplot(df_imp, aes(abs_loading)) +
    geom_histogram(bins = 30, color = "black") +
    geom_vline(xintercept = thr, linetype = 2) +
    theme_minimal(11) +
    labs(title = paste0("Umbral ", q, " (", name, ", TRAIN)"),
         x = "|loading| (max comp1–2 en TRAIN)", y = "Frecuencia")
  ggsave(file.path(dir_base, sprintf("hist_importancia_train_q%.2f.jpg", q)),
         p_hist_thr, width = 7, height = 5, dpi = 300)
  print(p_hist_thr)
  
  scores_tr <- as.data.frame(mdl2$variates$X)[, 1:2, drop = FALSE]
  colnames(scores_tr) <- c("comp1","comp2")
  df_tr <- cbind(scores_tr, grupo = as.factor(y_tr), id = rownames(scores_tr))
  p_scores_tr <- ggplot(df_tr, aes(comp1, comp2, color = grupo, label = id)) +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    theme_minimal(11) +
    labs(title = paste0("PLS-DA final (", name, ") · TRAIN"),
         x = paste0("PC1: ", var_pc[1], "%"),
         y = paste0("PC2: ", var_pc[2], "%"),
         color = "Grupo")+ scale_color_manual(values = c(FD="red", SO="blue", NP="black"))

  ggsave(file.path(dir_base, "plsda_scores_train.jpg"), p_scores_tr, width = 7, height = 5, dpi = 300)
  print(p_scores_tr)
  
  invisible(list(
    selected_vars = keep,
    threshold = thr,
    center = ctr_new, scale = scl_new,
    # -> orientación vars×muestras
    X_vs_train = X_vs_train,
    X_vs_test  = X_vs_test,
    # -> orientación muestras×vars  **NUEVO**
    X_sv_train = X_sv_train,
    X_sv_test  = X_sv_test,
    model_train_final = mdl2,
    train_ids = tr_ids,
    test_ids  = te_ids,
    paths = list(
      hist_threshold_train = file.path(dir_base, sprintf("hist_importancia_train_q%.2f.jpg", q)),
      scores_train         = file.path(dir_base, "plsda_scores_train.jpg")
    )
  ))
}

# === Extractor estándar para modelado: devuelve siempre muestras×vars alineado ===
extract_for_model <- function(fx_obj) {
  if (!all(c("X_sv_train","X_sv_test") %in% names(fx_obj))) {
    stop("No se hallaron matrices finales en fx_obj. Falta X_sv_train/X_sv_test.")
  }
  Xtr <- fx_obj$X_sv_train
  Xte <- fx_obj$X_sv_test
  if (!identical(colnames(Xtr), colnames(Xte)))
    stop("Columnas de train/test no coinciden.")
  list(
    X_train = Xtr,
    X_test  = Xte,
    tX_train = t(Xtr),
    tX_test  = t(Xte),
    center = fx_obj$center,
    scale  = fx_obj$scale,
    vars   = colnames(Xtr),
    train_ids = rownames(Xtr),
    test_ids  = rownames(Xte)
  )
}



##==codigo====

# 1) Metadata base
transcriptomica <- fun_process_transcriptomics("../data/transcriptomica.xlsx")

# 2) Split global único
sp <- define_split_ids(
  samples_metadata = transcriptomica$samples_metadata,
  grupo_col = "grupo",
  outlier_patterns = c("NP_02","NP_03","SO_14","6M"),
  p = 0.6, seed = 42
)

# 3) Proteoma (NO volver a quitar outliers aquí)
res_prot <- process_proteoma(
  path_citoquina = "../data/citokina.csv",
  path_miokina   = "../data/miokina.csv",
  metadata = transcriptomica$samples_metadata,
  y = sp$y_master,
  fixed_keep_ids  = sp$keep_ids,
  fixed_train_ids = sp$train_ids,
  fixed_test_ids  = sp$test_ids,
  remove_outliers = FALSE
)

# 4) Transcriptómica
res_tx <- prepare_transcriptomics_for_model(
  trans_obj = transcriptomica,
  metadata  = transcriptomica$samples_metadata,
  y = sp$y_master,
  fixed_keep_ids  = sp$keep_ids,
  fixed_train_ids = sp$train_ids,
  fixed_test_ids  = sp$test_ids,
  validate_with_bitr = TRUE,
  orgdb = "org.Hs.eg.db"
)

# 5) Metaboloma
res_met <- process_metaboloma(
  path_pos = "../data/5_03. METAHEALTH_POS_3 grupos_filtered_normalized.xlsx",
  path_neg = "../data/5_04. METAHEALTH_NEG_3 grupos_filtered_normalized.xlsx",
  metadata = transcriptomica$samples_metadata,
  y = sp$y_master,
  fixed_keep_ids  = sp$keep_ids,
  fixed_train_ids = sp$train_ids,
  fixed_test_ids  = sp$test_ids,
  align_by_suffix = TRUE,
  suffix_pattern  = ".*_(\\d+)"
)

# 6) Clínicos
res_cli <- process_clinicos(
  path_excel = "../data/clinicos.xlsx",
  metadata   = transcriptomica$samples_metadata,
  y = sp$y_master,
  fixed_keep_ids  = sp$keep_ids,
  fixed_train_ids = sp$train_ids,
  fixed_test_ids  = sp$test_ids,
  vars_clinicas   = CLINICAL_VARS
)

stopifnot(
  identical(rownames(res_prot$split$X_train_scaled), rownames(res_tx$split$X_train_scaled)),
  identical(rownames(res_prot$split$X_train_scaled), rownames(res_met$split$X_train_scaled)),
  identical(rownames(res_prot$split$X_train_scaled), rownames(res_cli$split$X_train_scaled)),
  identical(rownames(res_prot$split$X_test_scaled),  rownames(res_tx$split$X_test_scaled)),
  identical(rownames(res_prot$split$X_test_scaled),  rownames(res_met$split$X_test_scaled)),
  identical(rownames(res_prot$split$X_test_scaled),  rownames(res_cli$split$X_test_scaled))
)

# === Nuevos TEST (ejemplos) ===
# Proteoma:
# new_prot <- cbind(citoquina_new, miokina_new); rownames(new_prot) <- new_ids
# Xt_prot <- res_prot$transform_test(new_prot)$X_new_scaled
# Transcriptómica:
# new_tx <- t(new_data_by_gene_new)  # muestras x genes (Entrez como columnas)
# Xt_tx  <- res_tx$transform_test(new_tx)$X_new_scaled

# === Plots: usar shim para clínicos (desescalar X_log_all antes de pasarlo) ===
save_block_plots(res_tx,   out_dir = PLOTS_DIR, name = "Transcriptoma")
save_block_plots(res_prot, out_dir = PLOTS_DIR, name = "Proteoma")
save_block_plots(res_met,  out_dir = PLOTS_DIR, name = "Metaboloma")

save_block_plots(
  res_cli,
  out_dir = PLOTS_DIR,
  name = "Clinicos"
)
sapply(res_cli$split[c("X_train","X_test","X_train_scaled","X_test_scaled")], function(m){
  M <- as.matrix(m)
  c(nNA = sum(is.na(M)), nInf = sum(is.infinite(M)))
})

# === Importancia SOLO con TRAIN por ómica ===
# Mapeo Entrez -> símbolo para transcriptómica
id_to_label_tx <- {
  a <- transcriptomica$anno_by_entrez
  setNames(ifelse(is.na(a$GeneSymbol) | a$GeneSymbol == "", rownames(a), a$GeneSymbol),
           rownames(a))
}

tx_imp <- run_plsda_importance_train(
  res_tx,
  name = "Transcriptoma",
  out_dir = PLOTS_DIR,
  id_to_label = id_to_label_tx
)

pr_imp <- run_plsda_importance_train(
  res_prot,
  name = "Proteoma",
  out_dir = PLOTS_DIR
)

me_imp <- run_plsda_importance_train(
  res_met,
  name = "Metaboloma",
  out_dir = PLOTS_DIR
)

cl_imp <- run_plsda_importance_train(
  res_cli,
  name = "Clinicos",
  out_dir = PLOTS_DIR
)

imp_tx <- tx_imp$importance_train
imp_pr <- pr_imp$importance_train
imp_me <- me_imp$importance_train
imp_cl <- cl_imp$importance_train

# === Umbral por cuantil en TRAIN, refit y plots ===
fx_tx <- apply_importance_threshold_and_refit_from_train(
  res_tx,
  imp_tx,
  q = Q_TX,
  name = "Transcriptoma",
  out_dir = PLOTS_DIR
)

fx_pr <- apply_importance_threshold_and_refit_from_train(
  res_prot,
  imp_pr,
  q = Q_PR,
  name = "Proteoma",
  out_dir = PLOTS_DIR
)

fx_me <- apply_importance_threshold_and_refit_from_train(
  res_met,
  imp_me,
  q = Q_ME,
  name = "Metaboloma",
  out_dir = PLOTS_DIR
)

fx_cl <- apply_importance_threshold_and_refit_from_train(
  res_cli,
  imp_cl,
  q = Q_CL,
  name = "Clinicos",
  out_dir = PLOTS_DIR
)

# ============================================================
# EMPAQUETADO FINAL + AUDITORÍA FINAL
# Coherente con:
# - split global único entre ómicas
# - escalado aprendido solo en TRAIN
# - TEST escalado con parámetros de TRAIN
# - clínicos también filtrados por PLS-DA en TRAIN usando Q_CL
# ============================================================

# ------------------------------------------------------------
# 1) Extractor final para modelado / MOFA2
#    Entrada fx_*: muestras x variables
#    Salida lista: variables x sujetos
# ------------------------------------------------------------

extract_for_model <- function(fx_obj, view_name = "view", map_dir = OUTPUT_DIR) {
  
  if (!all(c("X_sv_train", "X_sv_test") %in% names(fx_obj))) {
    stop("No se hallaron matrices finales en fx_obj. Falta X_sv_train/X_sv_test.")
  }
  
  Xtr <- fx_obj$X_sv_train
  Xte <- fx_obj$X_sv_test
  
  Xtr <- as.matrix(Xtr)
  Xte <- as.matrix(Xte)
  mode(Xtr) <- "numeric"
  mode(Xte) <- "numeric"
  
  if (is.null(rownames(Xtr)) || is.null(rownames(Xte))) {
    stop("X_sv_train/X_sv_test deben tener rownames con IDs de muestra.")
  }
  
  if (is.null(colnames(Xtr)) || is.null(colnames(Xte))) {
    stop("X_sv_train/X_sv_test deben tener colnames con variables.")
  }
  
  stopifnot(identical(colnames(Xtr), colnames(Xte)))
  stopifnot(length(intersect(rownames(Xtr), rownames(Xte))) == 0)
  
  if (anyNA(Xtr) || anyNA(Xte)) {
    stop("Hay NA en matrices finales.")
  }
  
  if (any(!is.finite(Xtr)) || any(!is.finite(Xte))) {
    stop("Hay Inf/NaN en matrices finales.")
  }
  
  # ------------------------------------------------------------
  # Saneamiento de nombres de features para MOFA/modelado
  # ------------------------------------------------------------
  old_vars <- colnames(Xtr)
  new_vars <- ascii_feature_names(old_vars)
  
  feature_name_map <- data.frame(
    view = view_name,
    feature_original = old_vars,
    feature_ascii = new_vars,
    changed = old_vars != new_vars,
    stringsAsFactors = FALSE
  )
  
  if (!dir.exists(map_dir)) {
    dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  write.csv(
    feature_name_map,
    file.path(map_dir, paste0("feature_name_map_ascii_", view_name, ".csv")),
    row.names = FALSE
  )
  
  if (any(feature_name_map$changed)) {
    message(
      "Nombres saneados en vista ", view_name, ": ",
      sum(feature_name_map$changed), " modificados."
    )
  }
  
  colnames(Xtr) <- new_vars
  colnames(Xte) <- new_vars
  
  list(
    train = t(Xtr),          # variables x sujetos
    test  = t(Xte),          # variables x sujetos
    vars  = new_vars,
    ids_train = rownames(Xtr),
    ids_test  = rownames(Xte),
    feature_name_map = feature_name_map
  )
}

# ------------------------------------------------------------
# 2) Extraer matrices finales
# ------------------------------------------------------------

ex_tx <- extract_for_model(fx_tx, view_name = "tx")
ex_pr <- extract_for_model(fx_pr, view_name = "pr")
ex_me <- extract_for_model(fx_me, view_name = "me")
ex_cl <- extract_for_model(fx_cl, view_name = "cl")

# ------------------------------------------------------------
# 3) Verificación fuerte de IDs finales entre ómicas
# ------------------------------------------------------------

stopifnot(
  identical(ex_tx$ids_train, ex_pr$ids_train),
  identical(ex_tx$ids_train, ex_me$ids_train),
  identical(ex_tx$ids_train, ex_cl$ids_train),
  identical(ex_tx$ids_test,  ex_pr$ids_test),
  identical(ex_tx$ids_test,  ex_me$ids_test),
  identical(ex_tx$ids_test,  ex_cl$ids_test)
)

cat("\n============================================================\n")
cat("CHECK IDS FINALES ENTRE ÓMICAS\n")
cat("============================================================\n")
cat("Train IDs iguales entre ómicas: TRUE\n")
cat("Test IDs iguales entre ómicas : TRUE\n")
cat("n_train:", length(ex_cl$ids_train), "\n")
cat("n_test :", length(ex_cl$ids_test), "\n")
cat("Intersección train/test:", length(intersect(ex_cl$ids_train, ex_cl$ids_test)), "\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 4) Grupos train/test
# ------------------------------------------------------------

make_grupo_df <- function(obj) {
  
  extract_pref <- function(ids) sub("_.*", "", ids)
  
  list(
    train = data.frame(
      id = obj$ids_train,
      grupo = factor(
        extract_pref(obj$ids_train),
        levels = c("NP", "FD", "SO")
      )
    ),
    test = data.frame(
      id = obj$ids_test,
      grupo = factor(
        extract_pref(obj$ids_test),
        levels = c("NP", "FD", "SO")
      )
    )
  )
}

grupo_df <- make_grupo_df(ex_cl)

cat("\n============================================================\n")
cat("DISTRIBUCIÓN DE GRUPOS\n")
cat("============================================================\n")
cat("\nTRAIN:\n")
print(table(grupo_df$train$grupo, useNA = "ifany"))
cat("\nTEST:\n")
print(table(grupo_df$test$grupo, useNA = "ifany"))
cat("============================================================\n")

# ------------------------------------------------------------
# 5) Objeto final
# ------------------------------------------------------------

lista <- list(
  tx = ex_tx,
  pr = ex_pr,
  me = ex_me,
  cl = ex_cl,
  grupo_train = grupo_df$train,
  grupo_test  = grupo_df$test,
  features_metadata = transcriptomica$feature_metadata,
  feature_name_maps = list(
    tx = ex_tx$feature_name_map,
    pr = ex_pr$feature_name_map,
    me = ex_me$feature_name_map,
    cl = ex_cl$feature_name_map
  )
)

saveRDS(
  lista,
  file.path(OUTPUT_DIR, paste0("ready_for_modeling_", ANALYSIS_NAME, ".rds"))
)

run_config_preprocess <- data.frame(
  analysis_name = ANALYSIS_NAME,
  outdir_base = OUTDIR_BASE,
  output_dir = OUTPUT_DIR,
  plots_dir = PLOTS_DIR,
  
  q_tx = Q_TX,
  q_pr = Q_PR,
  q_me = Q_ME,
  q_cl = Q_CL,
  
  clinical_include_raw = CLINICAL_INCLUDE_RAW,
  clinical_exclude_raw = CLINICAL_EXCLUDE_RAW,
  clinical_vars_final = paste(CLINICAL_VARS, collapse = "|"),
  
  n_tx_selected = length(fx_tx$selected_vars),
  n_pr_selected = length(fx_pr$selected_vars),
  n_me_selected = length(fx_me$selected_vars),
  n_cl_selected = length(fx_cl$selected_vars),
  
  tx_selected = paste(fx_tx$selected_vars, collapse = "; "),
  pr_selected = paste(fx_pr$selected_vars, collapse = "; "),
  me_selected = paste(fx_me$selected_vars, collapse = "; "),
  cl_selected = paste(fx_cl$selected_vars, collapse = "; "),
  
  train_ids = paste(ex_cl$ids_train, collapse = "; "),
  test_ids  = paste(ex_cl$ids_test, collapse = "; "),
  
  stringsAsFactors = FALSE
)
write.csv(
  run_config_preprocess,
  file.path(OUTPUT_DIR, "run_config_preprocess.csv"),
  row.names = FALSE
)

saveRDS(
  run_config_preprocess,
  file.path(OUTPUT_DIR, "run_config_preprocess.rds")
)

cat("\n============================================================\n")
cat("CONFIGURACIÓN FINAL GUARDADA\n")
cat("============================================================\n")
print(run_config_preprocess)
cat("============================================================\n")
# ============================================================
# PLOTS PLS-DA SOBRE MATRICES FINALES
# ============================================================

group_col <- "grupo"
label_col <- "ID"
palette   <- c(FD = "red", SO = "blue", NP = "black")

grupo_vec <- setNames(
  as.factor(transcriptomica$samples_metadata$grupo),
  rownames(transcriptomica$samples_metadata)
)

plsda_plot <- function(X, meta_df, ttl,
                       group_col = "grupo",
                       label_col = NULL,
                       palette = c(FD = "red", SO = "blue", NP = "black")) {
  
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("mixOmics no instalado.")
  }
  
  y <- as.factor(meta_df[[group_col]])
  
  mdl <- mixOmics::plsda(X, y, ncomp = 2, scale = FALSE)
  var <- round(100 * mdl$prop_expl_var$X, 2)
  
  sc <- as.data.frame(mdl$variates$X)
  sc$grupo <- y
  
  if (!is.null(label_col) && label_col %in% colnames(meta_df)) {
    sc$label <- meta_df[[label_col]]
  } else {
    sc$label <- rownames(meta_df)
  }
  
  ggplot(sc, aes(comp1, comp2, color = grupo, label = label)) +
    geom_point(size = 2, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")
    ) +
    labs(
      title = ttl,
      x = paste0("PC1: ", var[1], "%"),
      y = paste0("PC2: ", var[2], "%"),
      color = "Grupo"
    ) +
    scale_color_manual(values = palette)
}

plot_plsda_split <- function(
    lst,
    split = c("train", "test"),
    out_dir = file.path(PLOTS_DIR, "PLSDA_final")
) {
  
  split <- match.arg(split)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (omic in c("tx", "pr", "me", "cl")) {
    
    M <- lst[[omic]][[split]]      # variables x sujetos
    ids <- colnames(M)             # sujetos
    
    meta_df <- data.frame(
      grupo = grupo_vec[ids],
      row.names = ids
    )
    
    X <- t(M)                      # muestras x variables
    
    ttl <- sprintf(
      "PLS-DA · %s · %s",
      toupper(omic),
      toupper(split)
    )
    
    p <- plsda_plot(X, meta_df, ttl)
    
    print(p)
    
    ggsave(
      filename = paste0("plsda_final_", omic, "_", split, ".jpg"),
      path = out_dir,
      plot = p,
      width = 7,
      height = 5,
      units = "in",
      dpi = 300,
      device = "jpeg",
      quality = 95
    )
  }
}
plot_plsda_split(lista, "train", out_dir = file.path(PLOTS_DIR, "PLSDA_final"))
plot_plsda_split(lista, "test",  out_dir = file.path(PLOTS_DIR, "PLSDA_final"))

# ============================================================
# AUDITORÍA FINAL DE ESCALADO PRE-MOFA2 / MODELADO
# No guarda nada. Solo imprime por consola.
# ============================================================

audit_scaled_matrix <- function(X, label = "matriz", tol_mean = 1e-10, tol_sd = 1e-10) {
  
  X <- as.matrix(X)
  mode(X) <- "numeric"
  
  n_na  <- sum(is.na(X))
  n_inf <- sum(is.infinite(X))
  n_nan <- sum(is.nan(X))
  
  means <- colMeans(X)
  sds   <- apply(X, 2, sd)
  
  tab <- data.frame(
    variable = colnames(X),
    mean = as.numeric(means),
    sd = as.numeric(sds),
    abs_mean = abs(as.numeric(means)),
    abs_sd_minus_1 = abs(as.numeric(sds) - 1),
    mean_ok = abs(as.numeric(means)) <= tol_mean,
    sd_ok = abs(as.numeric(sds) - 1) <= tol_sd,
    stringsAsFactors = FALSE
  )
  
  cat("\n============================================================\n")
  cat("AUDITORÍA:", label, "\n")
  cat("Dimensión:", nrow(X), "muestras x", ncol(X), "variables\n")
  cat("NA:", n_na, "| NaN:", n_nan, "| Inf:", n_inf, "\n")
  cat("Máx |media|:", max(tab$abs_mean, na.rm = TRUE), "\n")
  cat("Máx |SD - 1|:", max(tab$abs_sd_minus_1, na.rm = TRUE), "\n")
  cat("Variables con media OK:", sum(tab$mean_ok), "/", nrow(tab), "\n")
  cat("Variables con SD OK:", sum(tab$sd_ok), "/", nrow(tab), "\n")
  cat("============================================================\n")
  
  print(tab[order(-tab$abs_mean, -tab$abs_sd_minus_1), ], row.names = FALSE)
  
  invisible(tab)
}

audit_test_scaled_with_train_params <- function(X_test_scaled, label = "TEST") {
  
  X <- as.matrix(X_test_scaled)
  mode(X) <- "numeric"
  
  tab <- data.frame(
    variable = colnames(X),
    mean_test = as.numeric(colMeans(X)),
    sd_test = as.numeric(apply(X, 2, sd)),
    min_test = as.numeric(apply(X, 2, min)),
    max_test = as.numeric(apply(X, 2, max)),
    stringsAsFactors = FALSE
  )
  
  cat("\n============================================================\n")
  cat("AUDITORÍA:", label, "\n")
  cat("Nota: TEST fue escalado con center/scale de TRAIN.\n")
  cat("Por eso TEST NO tiene obligación de tener media 0 ni SD 1.\n")
  cat("Dimensión:", nrow(X), "muestras x", ncol(X), "variables\n")
  cat("NA:", sum(is.na(X)), "| NaN:", sum(is.nan(X)), "| Inf:", sum(is.infinite(X)), "\n")
  cat("============================================================\n")
  
  print(tab, row.names = FALSE)
  
  invisible(tab)
}

audit_split_ids <- function(res_list) {
  
  cat("\n============================================================\n")
  cat("AUDITORÍA DE IDS TRAIN/TEST ENTRE ÓMICAS\n")
  cat("============================================================\n")
  
  train_ids <- lapply(res_list, function(x) rownames(x$split$X_train_scaled))
  test_ids  <- lapply(res_list, function(x) rownames(x$split$X_test_scaled))
  
  ref_train <- train_ids[[1]]
  ref_test  <- test_ids[[1]]
  
  for (nm in names(res_list)) {
    cat("\n[", nm, "]\n", sep = "")
    cat("n_train:", length(train_ids[[nm]]), "\n")
    cat("n_test :", length(test_ids[[nm]]), "\n")
    cat("Train igual al primer bloque:", identical(ref_train, train_ids[[nm]]), "\n")
    cat("Test igual al primer bloque :", identical(ref_test, test_ids[[nm]]), "\n")
    cat("Intersección train/test:", length(intersect(train_ids[[nm]], test_ids[[nm]])), "\n")
  }
  
  invisible(list(train_ids = train_ids, test_ids = test_ids))
}

audit_fx_object <- function(fx, label = "fx") {
  
  cat("\n============================================================\n")
  cat("AUDITORÍA OBJETO FINAL:", label, "\n")
  cat("============================================================\n")
  
  Xtr <- fx$X_sv_train
  Xte <- fx$X_sv_test
  
  if (is.null(Xtr) || is.null(Xte)) {
    stop(label, " no contiene X_sv_train / X_sv_test.")
  }
  
  cat("Dim train:", nrow(Xtr), "x", ncol(Xtr), "\n")
  cat("Dim test :", nrow(Xte), "x", ncol(Xte), "\n")
  cat("Columnas train/test idénticas:", identical(colnames(Xtr), colnames(Xte)), "\n")
  cat("Intersección IDs train/test:", length(intersect(rownames(Xtr), rownames(Xte))), "\n")
  cat("N variables seleccionadas:", ncol(Xtr), "\n")
  
  audit_scaled_matrix(Xtr, paste0(label, " · TRAIN final"))
  audit_test_scaled_with_train_params(Xte, paste0(label, " · TEST final"))
  
  invisible(TRUE)
}

# ------------------------------------------------------------
# 6) Auditoría de splits originales
# ------------------------------------------------------------

res_omics <- list(
  Transcriptoma = res_tx,
  Proteoma      = res_prot,
  Metaboloma    = res_met,
  Clinicos      = res_cli
)

audit_ids <- audit_split_ids(res_omics)

# ------------------------------------------------------------
# 7) Auditoría de matrices escaladas originales
# ------------------------------------------------------------

audit_original <- list()

for (nm in names(res_omics)) {
  
  res_now <- res_omics[[nm]]
  
  audit_original[[nm]] <- list(
    train = audit_scaled_matrix(
      res_now$split$X_train_scaled,
      paste0(nm, " · TRAIN scaled")
    ),
    test = audit_test_scaled_with_train_params(
      res_now$split$X_test_scaled,
      paste0(nm, " · TEST scaled con parámetros de TRAIN")
    )
  )
}

# ------------------------------------------------------------
# 8) Auditoría de objetos finales fx_*
# ------------------------------------------------------------

fx_list <- list(
  Transcriptoma = fx_tx,
  Proteoma      = fx_pr,
  Metaboloma    = fx_me,
  Clinicos      = fx_cl
)

for (nm in names(fx_list)) {
  audit_fx_object(fx_list[[nm]], nm)
}

# ------------------------------------------------------------
# 9) Resumen compacto final
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("RESUMEN FINAL DE AUDITORÍA\n")
cat("============================================================\n")

for (nm in names(res_omics)) {
  
  Xtr <- as.matrix(res_omics[[nm]]$split$X_train_scaled)
  Xte <- as.matrix(res_omics[[nm]]$split$X_test_scaled)
  
  cat("\n[", nm, "]\n", sep = "")
  cat("TRAIN max |media|     :", max(abs(colMeans(Xtr))), "\n")
  cat("TRAIN max |SD - 1|    :", max(abs(apply(Xtr, 2, sd) - 1)), "\n")
  cat("TRAIN NA/Inf          :", sum(is.na(Xtr)), "/", sum(is.infinite(Xtr)), "\n")
  cat("TEST  NA/Inf          :", sum(is.na(Xte)), "/", sum(is.infinite(Xte)), "\n")
  cat("Train/Test disjoint   :", length(intersect(rownames(Xtr), rownames(Xte))) == 0, "\n")
}

cat("\n============================================================\n")
cat("INTERPRETACIÓN\n")
cat("============================================================\n")
cat("Si TRAIN tiene media ~0 y SD ~1, el escalado fue aprendido correctamente.\n")
cat("Si TEST no tiene media 0 ni SD 1, eso es normal: TEST se escaló con parámetros de TRAIN.\n")
cat("Si no hay NA/Inf y los IDs train/test son disjuntos e iguales entre ómicas, la entrada está coherente.\n")
cat("============================================================\n")
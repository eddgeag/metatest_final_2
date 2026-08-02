library(MOFA2)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(pheatmap)
library(factoextra)
library(tibble)
library(ggrepel)
library(Hmisc)   
library(limma)
library(Bolstad)
library(tidyr)
library(ggnewscale)


## ============================================================
## 0) CONFIGURACIÓN PARAMETRIZADA
## Cambia SOLO ANALYSIS_NAME y PANEL_CHOICE para cada corrida
## ============================================================
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

ANALYSIS_NAME <- get_env_chr("ANALYSIS_NAME", required = TRUE)
OUTDIR_BASE   <- get_env_chr("OUTDIR_BASE", ".")

## ============================================================
## Qué diagnóstico ejecutar
## BAYES_INPUT = diagnostica las features finales del modelamiento BRMS/MOFA
## PANEL_FINAL = diagnostica el panel final tras multimodelo
## BOTH        = ejecuta ambos diagnósticos
## ============================================================

DIAG_TARGET <- toupper(get_env_chr("DIAG_TARGET", "BOTH"))

DIAG_TARGET <- dplyr::recode(
  DIAG_TARGET,
  "BAYES"        = "BAYES_INPUT",
  "BRMS"         = "BAYES_INPUT",
  "BRMS_INPUT"   = "BAYES_INPUT",
  "MOFA_BRMS"    = "BAYES_INPUT",
  "PANEL"        = "PANEL_FINAL",
  "FINAL"        = "PANEL_FINAL",
  .default = DIAG_TARGET
)

if (!DIAG_TARGET %in% c("BAYES_INPUT", "PANEL_FINAL", "BOTH")) {
  stop(
    "DIAG_TARGET debe ser BAYES_INPUT, PANEL_FINAL o BOTH. Valor actual: ",
    DIAG_TARGET
  )
}

## Solo aplica cuando DIAG_TARGET = PANEL_FINAL
## FINAL = panel final MOFA + >=2/3 modelos supervisados
## ALL   = tabla completa de consenso multimodelo
PANEL_CHOICE <- toupper(get_env_chr("PANEL_CHOICE", "FINAL"))

PANEL_CHOICE <- dplyr::recode(
  PANEL_CHOICE,
  "RELAXED"   = "FINAL",
  "STRICT"    = "FINAL",
  "M2OF3"     = "FINAL",
  "MOFA_2OF3" = "FINAL",
  .default = PANEL_CHOICE
)

if (!PANEL_CHOICE %in% c("FINAL", "ALL")) {
  stop(
    "PANEL_CHOICE debe ser FINAL o ALL. Valor actual: ",
    PANEL_CHOICE
  )
}
## Qué modelamiento usar como entrada
## Opciones:
##   NEW  = modelo actual: ANALYSIS_DIR/modelamiento/
##   OLD  = modelo anterior: ANALYSIS_DIR/modelado_anterior/
##   BOTH = relanza este mismo script para NEW y OLD
MODEL_SOURCE <- toupper(get_env_chr("MODEL_SOURCE", "NEW"))

if (!MODEL_SOURCE %in% c("NEW", "OLD", "BOTH")) {
  stop(
    "MODEL_SOURCE debe ser NEW, OLD o BOTH. Valor actual: ",
    MODEL_SOURCE
  )
}
ANALYSIS_DIR <- file.path(OUTDIR_BASE, ANALYSIS_NAME)

## Entradas generadas por scripts anteriores
RFM_FILE <- file.path(
  ANALYSIS_DIR,
  paste0("ready_for_modeling_", ANALYSIS_NAME, ".rds")
)

MOFA_MODEL_FILE <- file.path(
  ANALYSIS_DIR,
  "integracion",
  "rds",
  paste0("modelo_mofa_final_", ANALYSIS_NAME, ".rds")
)

## ============================================================
## Selección de entrada downstream: modelo nuevo vs modelo anterior
## ============================================================

## ============================================================
## Relanzar combinaciones:
##   MODEL_SOURCE = NEW / OLD
##   DIAG_TARGET  = BAYES_INPUT / PANEL_FINAL
## ============================================================

if (MODEL_SOURCE == "BOTH" || DIAG_TARGET == "BOTH") {
  
  cmd_args <- commandArgs(FALSE)
  script_file <- sub("^--file=", "", grep("^--file=", cmd_args, value = TRUE)[1])
  
  if (is.na(script_file) || !nzchar(script_file)) {
    stop(
      "MODEL_SOURCE=BOTH o DIAG_TARGET=BOTH requiere ejecutar con Rscript script.R. ",
      "Si usas source(), corre cada combinación por separado."
    )
  }
  
  src_vec <- if (MODEL_SOURCE == "BOTH") c("NEW", "OLD") else MODEL_SOURCE
  tgt_vec <- if (DIAG_TARGET == "BOTH") c("BAYES_INPUT", "PANEL_FINAL") else DIAG_TARGET
  
  for (src in src_vec) {
    for (tgt in tgt_vec) {
      
      cat("\n============================================================\n")
      cat("RELANZANDO DOWNSTREAM | MODEL_SOURCE =", src, "| DIAG_TARGET =", tgt, "\n")
      cat("============================================================\n")
      
      status <- system2(
        command = file.path(R.home("bin"), "Rscript"),
        args = shQuote(script_file),
        env = c(
          paste0("ANALYSIS_NAME=", ANALYSIS_NAME),
          paste0("OUTDIR_BASE=", OUTDIR_BASE),
          paste0("MODEL_SOURCE=", src),
          paste0("DIAG_TARGET=", tgt),
          paste0("PANEL_CHOICE=", PANEL_CHOICE)
        )
      )
      
      if (!identical(status, 0L)) {
        stop(
          "Falló downstream para MODEL_SOURCE = ",
          src,
          " | DIAG_TARGET = ",
          tgt
        )
      }
    }
  }
  
  quit(save = "no", status = 0)
}

MODEL_DIR_NAME <- dplyr::case_when(
  MODEL_SOURCE == "NEW" ~ "modelamiento",
  MODEL_SOURCE == "OLD" ~ "modelado_anterior",
  TRUE ~ NA_character_
)

PANEL_DIR_NAME <- dplyr::case_when(
  MODEL_SOURCE == "NEW" ~ "panel_final_multimodelo",
  MODEL_SOURCE == "OLD" ~ "panel_final_multimodelo_anterior",
  TRUE ~ NA_character_
)


DIAG_DIR_NAME <- dplyr::case_when(
  MODEL_SOURCE == "NEW" & DIAG_TARGET == "BAYES_INPUT" ~
    "diagnostico_bayes_features",
  
  MODEL_SOURCE == "OLD" & DIAG_TARGET == "BAYES_INPUT" ~
    "diagnostico_modelado_anterior_bayes_features",
  
  MODEL_SOURCE == "NEW" & DIAG_TARGET == "PANEL_FINAL" ~
    paste0("diagnostico_panel_", tolower(PANEL_CHOICE)),
  
  MODEL_SOURCE == "OLD" & DIAG_TARGET == "PANEL_FINAL" ~
    paste0("diagnostico_modelado_anterior_panel_", tolower(PANEL_CHOICE)),
  
  TRUE ~ NA_character_
)
MODELAMIENTO_RDS_DIR <- file.path(
  ANALYSIS_DIR,
  MODEL_DIR_NAME,
  "rds"
)

MODELAMIENTO_TABLES_DIR <- file.path(
  ANALYSIS_DIR,
  MODEL_DIR_NAME,
  "tables"
)

PANEL_TABLES_DIR <- file.path(
  ANALYSIS_DIR,
  PANEL_DIR_NAME,
  "tables"
)
PANEL_RDS_DIR <- file.path(
  ANALYSIS_DIR,
  PANEL_DIR_NAME,
  "rds"
)

PANEL_OBJECT_FILE <- file.path(
  PANEL_RDS_DIR,
  "panel_final_multimodel_object.rds"
)

## Salidas del modelamiento bayesiano inicial.
## Estas son la entrada del multimodelo.
BAYES_FEATURES_RDS_FILE <- file.path(
  MODELAMIENTO_RDS_DIR,
  "features_finales_integradas_mofa_brms.rds"
)

BAYES_FEATURES_CSV_FILE <- file.path(
  MODELAMIENTO_TABLES_DIR,
  "features_finales_integradas_mofa_brms.csv"
)

BAYES_FEATURES_OBJECT_FILE <- file.path(
  MODELAMIENTO_RDS_DIR,
  "features_finales_object.rds"
)


## Salidas de este script
DIAG_DIR <- file.path(
  ANALYSIS_DIR,
  DIAG_DIR_NAME
)
DIAG_TABLES_DIR <- file.path(DIAG_DIR, "tables")
DIAG_RDS_DIR    <- file.path(DIAG_DIR, "rds")

dir.create(DIAG_DIR,        recursive = TRUE, showWarnings = FALSE)
dir.create(DIAG_TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAG_RDS_DIR,    recursive = TRUE, showWarnings = FALSE)

out_diag <- function(...) {
  file.path(DIAG_DIR, ...)
}

out_table <- function(...) {
  file.path(DIAG_TABLES_DIR, ...)
}

out_rds <- function(...) {
  file.path(DIAG_RDS_DIR, ...)
}

dir_create <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

cat("\n============================================================\n")
cat("DIAGNÓSTICO PARAMETRIZADO\n")
cat("============================================================\n")
cat("ANALYSIS_NAME:", ANALYSIS_NAME, "\n")
cat("OUTDIR_BASE:", OUTDIR_BASE, "\n")
cat("PANEL_CHOICE:", PANEL_CHOICE, "\n")
cat("MODEL_SOURCE:", MODEL_SOURCE, "\n")
cat("MODEL_DIR_NAME:", MODEL_DIR_NAME, "\n")
cat("PANEL_DIR_NAME:", PANEL_DIR_NAME, "\n")
cat("ANALYSIS_DIR:", ANALYSIS_DIR, "\n")
cat("RFM_FILE:", RFM_FILE, "\n")
cat("MOFA_MODEL_FILE:", MOFA_MODEL_FILE, "\n")
cat("MODELAMIENTO_RDS_DIR:", MODELAMIENTO_RDS_DIR, "\n")
cat("PANEL_TABLES_DIR:", PANEL_TABLES_DIR, "\n")
cat("DIAG_DIR:", DIAG_DIR, "\n")
cat("DIAG_TARGET:", DIAG_TARGET, "\n")
cat("PANEL_RDS_DIR:", PANEL_RDS_DIR, "\n")
cat("PANEL_OBJECT_FILE:", PANEL_OBJECT_FILE, "\n")
cat("BAYES_FEATURES_RDS_FILE:", BAYES_FEATURES_RDS_FILE, "\n")
cat("BAYES_FEATURES_CSV_FILE:", BAYES_FEATURES_CSV_FILE, "\n")
cat("BAYES_FEATURES_OBJECT_FILE:", BAYES_FEATURES_OBJECT_FILE, "\n")
cat("============================================================\n")


plot_pca_feats <- function(pcx, df_weights_unique, feats, axes = c(1,2)) {
  
  # Paleta fija
  omics_colors <- c(
    "Transcriptomics" = "red",
    "Proteomics"      = "black",
    "Metabolomics"    = "olivedrab",
    "Clinical"        = "blue"
  )
  
  # Extraer loadings de PCA
  loadings <- as.data.frame(pcx$rotation[, axes]) %>%
    tibble::rownames_to_column("Feature") %>%
    dplyr::filter(Feature %in% feats) %>%
    dplyr::left_join(df_weights_unique, by = "Feature") %>%
    dplyr::filter(!is.na(view)) %>%
    dplyr::mutate(
      label_pretty = dplyr::coalesce(label_pretty, Feature)
    )
  # Renombrar ejes
  axis_names <- paste0("PC", axes)
  colnames(loadings)[2:3] <- axis_names
  
  varianza <-  round(100*(pcx$sdev^2)/sum(pcx$sdev^2),2)
  # --- Plot ---
  ggplot(loadings, aes_string(
    x = axis_names[1],
    y = axis_names[2],
    label = "label_pretty",
    color = "view")) +
    geom_point(size = 3) +
    geom_text_repel(max.overlaps = 50, size = 3) +
    scale_color_manual(values = omics_colors) +
    labs(
      title = "PCA loadings",
      x = paste0(axis_names[1], " (", varianza[as.numeric(gsub("PC","",axis_names[1]))], "%)"),
      y = paste0(axis_names[2], " (", varianza[as.numeric(gsub("PC","",axis_names[2]))], "%)"),
      color = "Ómic"
    ) +
    theme_minimal(base_size = 14)
}


get_corr <- function(df, feats) {
  
  feats <- intersect(feats, colnames(df))
  
  if (!all(c("Factor1", "Factor2") %in% colnames(df))) {
    stop("Faltan Factor1 y/o Factor2 en el dataframe.")
  }
  
  if (length(feats) == 0) {
    stop("No hay features válidas para correlacionar con Factor1/Factor2.")
  }
  
  stats::cor(
    df[, c("Factor1", "Factor2"), drop = FALSE],
    df[, feats, drop = FALSE],
    method = "spearman",
    use = "pairwise.complete.obs"
  )
}
make_corr_heatmap_global <- function(data, annotation_col, ann_colors) {
  
  sub <- data %>%
    dplyr::select(-dplyr::any_of("y"))
  
  feats <- intersect(annotation_col$Label, colnames(sub))
  
  if (!all(c("Factor1", "Factor2") %in% colnames(sub))) {
    stop("Faltan Factor1 y/o Factor2 en all_data para correlación global.")
  }
  
  if (length(feats) == 0) {
    stop("No hay features comunes entre annotation_col$Label y colnames(all_data).")
  }
  
  ann_sub <- annotation_col[feats, , drop = FALSE]
  ann_sub <- ann_sub[order(ann_sub$View), , drop = FALSE]
  feats <- rownames(ann_sub)
  
  corr_sub <- get_corr(sub, feats)
  colnames(corr_sub) <- feats
  pretty_labels <- annotation_col[feats, "Label_pretty"]
  pretty_labels <- as.character(pretty_labels)
  pretty_labels[is.na(pretty_labels)] <- feats
  pheatmap::pheatmap(
    corr_sub,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 12,
    fontsize_col = 7,
    main = "Global Spearman correlation (Factors vs Features)",
    annotation_col = ann_sub[, "View", drop = FALSE],
    annotation_colors = ann_colors,
    border_color = NA,
    labels_col = pretty_labels
  )
}

make_corr_heatmap_group <- function(data,
                                    group_name,
                                    annotation_col,
                                    ann_colors) {
  
  sub <- data %>%
    dplyr::filter(.data$y == group_name) %>%
    dplyr::select(-dplyr::any_of("y"))
  
  feats <- intersect(annotation_col$Label, colnames(sub))
  
  if (!all(c("Factor1", "Factor2") %in% colnames(sub))) {
    stop("Faltan Factor1 y/o Factor2 para el grupo: ", group_name)
  }
  
  if (length(feats) == 0) {
    stop("No hay features comunes para el grupo: ", group_name)
  }
  
  ann_sub <- annotation_col[feats, , drop = FALSE]
  ann_sub <- ann_sub[order(ann_sub$View), , drop = FALSE]
  feats <- rownames(ann_sub)
  
  corr_sub <- get_corr(sub, feats)
  colnames(corr_sub) <- feats
  pretty_labels <- annotation_col[feats, "Label_pretty"]
  pretty_labels <- as.character(pretty_labels)
  pretty_labels[is.na(pretty_labels)] <- feats
  pheatmap::pheatmap(
    corr_sub,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 12,
    fontsize_col = 7,
    main = paste("Spearman correlation:", group_name),
    annotation_col = ann_sub[, "View", drop = FALSE],
    annotation_colors = ann_colors,
    labels_col = pretty_labels,
    border_color = NA
  )
}


normalize_name <- function(x) {
  
  x %>%
    gsub("_", "", .) %>%
    gsub("-", "", .) %>%
    gsub("\\.", "", .) %>%
    toupper()
}

get_corr_sig <- function(df, feats) {
  mat <- as.matrix(df[, feats, drop = FALSE])
  res <- Hmisc::rcorr(mat, type = "spearman")
  corr <- res$r
  pval <- res$P
  sig <- ifelse(pval < 0.05, "*", "")
  diag(corr) <- 1
  diag(sig)  <- ""
  list(corr = corr, sig = sig)
  
}

# --- Heatmap sin dendrograma ---
make_feature_corr_heatmap <- function(data, annotation_col, ann_colors, group_name = "Global") {
  
  if (group_name != "Global") {
    sub <- data %>%
      dplyr::filter(.data$y == group_name)
  } else {
    sub <- data
  }
  
  drop_cols <- c("y", "Factor1", "Factor2")
  
  sub <- sub %>%
    dplyr::select(-dplyr::any_of(drop_cols))
  
  feats <- intersect(annotation_col$Label, colnames(sub))
  
  if (length(feats) < 2) {
    stop(
      "Muy pocas features comunes para heatmap de correlación en grupo: ",
      group_name,
      ". Features comunes = ",
      length(feats)
    )
  }
  
  ann_sub <- annotation_col[feats, , drop = FALSE]
  ann_sub <- ann_sub[order(ann_sub$View), , drop = FALSE]
  feats <- rownames(ann_sub)
  
  corr_res <- get_corr_sig(sub, feats)
  pretty_labels <- annotation_col[feats, "Label_pretty"]
  pretty_labels <- as.character(pretty_labels)
  pretty_labels[is.na(pretty_labels)] <- feats
  names(pretty_labels) <- feats
  pheatmap::pheatmap(
    corr_res$corr[feats, feats],
    color = colorRampPalette(c("blue","white","red"))(100),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    main = paste("Spearman correlations between panel features -", group_name),
    annotation_col = ann_sub[, "View", drop = FALSE],
    annotation_colors = ann_colors,
    labels_row = pretty_labels,
    labels_col = pretty_labels,
    display_numbers = corr_res$sig[feats, feats],
    number_color = "black",
    border_color = NA
  )
}

contrast_t_test <- function(datas, factor_, response, contrast_expr) {
  
  if (!all(c(factor_, response) %in% colnames(datas))) {
    stop("El factor y la variable de respuesta deben estar en las columnas del dataframe.")
  }
  
  datas <- datas %>%
    dplyr::select(dplyr::all_of(c(factor_, response))) %>%
    dplyr::filter(!is.na(.data[[factor_]]), !is.na(.data[[response]]))
  
  datas[[factor_]] <- factor(datas[[factor_]])
  
  contrast_clean <- gsub(" ", "", contrast_expr)
  groups <- unlist(strsplit(contrast_clean, "-"))
  
  if (length(groups) != 2) {
    warning("Contraste no reconocido: ", contrast_expr)
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  
  g1 <- groups[1]
  g2 <- groups[2]
  
  if (!all(c(g1, g2) %in% levels(datas[[factor_]]))) {
    warning("Los grupos del contraste no están en el factor: ", contrast_expr)
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  
  x <- datas[[response]][datas[[factor_]] == g1]
  y <- datas[[response]][datas[[factor_]] == g2]
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) < 2 || length(y) < 2) {
    warning("Muy pocas observaciones para ", response, " en contraste ", contrast_expr)
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  
  vx <- stats::var(x)
  vy <- stats::var(y)
  
  if (!is.finite(vx) || !is.finite(vy) || vx == 0 || vy == 0) {
    warning("Varianza cero/no finita para ", response, " en contraste ", contrast_expr)
    
    tt <- tryCatch(
      stats::t.test(x, y, var.equal = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(tt)) {
      return(list(statistic = NA_real_, p.value = NA_real_))
    }
    
    return(list(
      statistic = unname(tt$statistic),
      p.value = unname(tt$p.value)
    ))
  }
  
  m_x <- mean(x)
  m_y <- mean(y)
  
  m <- c(m_x, m_y)
  n0 <- c(1 / vx, 1 / vy)
  sig.med <- median(c(sd(x), sd(y)))
  kappa <- 1
  
  res <- tryCatch(
    Bolstad::bayes.t.test(
      x,
      y,
      var.equal = FALSE,
      prior = "joint.conj",
      m = m,
      n0 = n0,
      sig.med = sig.med,
      kappa = kappa
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    return(list(
      statistic = unname(res$statistic),
      p.value = unname(res$p.value)
    ))
  }
  
  ## Fallback Welch si bayes.t.test falla
  tt <- tryCatch(
    stats::t.test(x, y, var.equal = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(tt)) {
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  
  list(
    statistic = unname(tt$statistic),
    p.value = unname(tt$p.value)
  )
}


plot_heatmap_stats_ordered <- function(statistic, p_adj, df_weights, alpha = 0.05) {
  
  # df_weights <- df_weights_unique
  # alpha <-  0.5
  # Pasar a formato largo
  stat_long <- statistic %>%
    as.data.frame() %>%
    rownames_to_column("Feature") %>%
    pivot_longer(-Feature, names_to = "Contrast", values_to = "Statistic")
  
  pval_long <- p_adj %>%
    as.data.frame() %>%
    rownames_to_column("Feature") %>%
    pivot_longer(-Feature, names_to = "Contrast", values_to = "p_adj")
  
  df_long <- left_join(stat_long, pval_long, by = c("Feature", "Contrast"))
  df_long <- inner_join(df_long, df_weights, by = "Feature")
  df_long <- df_long %>%
    dplyr::mutate(
      sig = dplyr::case_when(
        is.na(p_adj) ~ "",
        p_adj < alpha ~ "*",
        TRUE ~ ""
      ),
      view = factor(
        view,
        levels = c("Transcriptomics", "Proteomics", "Metabolomics", "Clinical")
      ),
      Feature = factor(Feature, levels = df_weights$Feature),
      label_pretty = dplyr::coalesce(label_pretty, as.character(Feature))
    )
  
  label_map <- df_weights %>%
    dplyr::distinct(Feature, label_pretty) %>%
    tibble::deframe()
  # Colores de ómicas
  omics_colors <- c(
    "Transcriptomics" = "red",
    "Proteomics"      = "black",
    "Metabolomics"    = "olivedrab",
    "Clinical"        = "blue"
  )
  
  # Data frame auxiliar para anotación de ómica
  df_bar <- df_long %>% distinct(Feature, view)
  
  ggplot(df_long, aes(x = Contrast, y = Feature)) +
    # Heatmap principal
    geom_tile(aes(fill = Statistic), color = "white") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = "Statistic"
    ) +
    
    # Estrellitas
    geom_text(aes(label = sig), color = "black", size = 3) +
    
    # Barra lateral con colores de ómica (otra escala de fill, independiente)
    new_scale_fill() +  # <- requiere library(ggnewscale)
    geom_tile(
      data = df_bar,
      aes(x = 0, y = Feature, fill = view),
      inherit.aes = FALSE,
      width = 0.3
    ) +
    scale_fill_manual(values = omics_colors, name = "Omic") +
    scale_y_discrete(labels = label_map) +
    facet_grid(view ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10,colour = "black",face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text.y = element_text(angle = 0, face = "bold",colour = "black"),
      panel.spacing.y = unit(0.5, "lines")
    ) +
    labs(
      title = "Empirical Bayes t-test",
      x = "Contrast", y = "Feature"
    )
}

## ============================================================
## CARGA DE OBJETOS DEL ANÁLISIS
## ============================================================

if (!file.exists(RFM_FILE)) {
  stop("No existe RFM_FILE: ", RFM_FILE)
}

if (!file.exists(MOFA_MODEL_FILE)) {
  stop("No existe MOFA_MODEL_FILE: ", MOFA_MODEL_FILE)
}

DF_TRAIN_FILE <- file.path(MODELAMIENTO_RDS_DIR, "df_train_final.rds")
DF_TEST_FILE  <- file.path(MODELAMIENTO_RDS_DIR, "df_test_final.rds")

if (!file.exists(DF_TRAIN_FILE)) {
  stop("No existe df_train_final.rds en: ", DF_TRAIN_FILE)
}

if (!file.exists(DF_TEST_FILE)) {
  stop("No existe df_test_final.rds en: ", DF_TEST_FILE)
}

## ============================================================
## Cargar salida del modelamiento del panel
## Preferencia:
##   1) panel_final_multimodel_object.rds
##   2) CSVs de tables/
## ============================================================

## ============================================================
## Cargar set de features para diagnóstico
##   BAYES_INPUT = features finales del modelamiento BRMS/MOFA
##   PANEL_FINAL = panel final del modelamiento multimodelo
## ============================================================

load_diagnostic_feature_set <- function(diag_target, panel_choice = "FINAL") {
  
  ## ------------------------------------------------------------
  ## 1) Diagnóstico de entrada: features finales BRMS/MOFA
  ## ------------------------------------------------------------
  if (diag_target == "BAYES_INPUT") {
    
    if (file.exists(BAYES_FEATURES_RDS_FILE)) {
      
      cat("\nFeatures bayesianas cargadas desde RDS:\n")
      cat(BAYES_FEATURES_RDS_FILE, "\n")
      
      feature_tbl <- readRDS(BAYES_FEATURES_RDS_FILE) %>%
        tibble::as_tibble()
      
      return(list(
        feature_tbl = feature_tbl,
        feature_file_used = BAYES_FEATURES_RDS_FILE,
        feature_input_type = "bayes_features_rds",
        diagnostic_source = "BAYES_INPUT"
      ))
    }
    
    if (file.exists(BAYES_FEATURES_CSV_FILE)) {
      
      cat("\nFeatures bayesianas cargadas desde CSV:\n")
      cat(BAYES_FEATURES_CSV_FILE, "\n")
      
      feature_tbl <- read.csv(BAYES_FEATURES_CSV_FILE, check.names = FALSE) %>%
        tibble::as_tibble()
      
      return(list(
        feature_tbl = feature_tbl,
        feature_file_used = BAYES_FEATURES_CSV_FILE,
        feature_input_type = "bayes_features_csv",
        diagnostic_source = "BAYES_INPUT"
      ))
    }
    
    stop(
      "No encuentro features finales del modelamiento bayesiano:\n",
      BAYES_FEATURES_RDS_FILE, "\n",
      BAYES_FEATURES_CSV_FILE, "\n",
      "Primero corre el script de modelamiento bayesiano."
    )
  }
  
  ## ------------------------------------------------------------
  ## 2) Diagnóstico de salida: panel final multimodelo
  ## ------------------------------------------------------------
  if (diag_target == "PANEL_FINAL") {
    
    if (file.exists(PANEL_OBJECT_FILE)) {
      
      cat("\nPanel cargado desde objeto RDS del modelamiento multimodelo:\n")
      cat(PANEL_OBJECT_FILE, "\n")
      
      panel_obj <- readRDS(PANEL_OBJECT_FILE)
      
      feature_tbl <- switch(
        panel_choice,
        "FINAL" = panel_obj$panel_final,
        "ALL"   = panel_obj$panel_feature_tbl,
        stop("PANEL_CHOICE no reconocido: ", panel_choice)
      )
      
      if (is.null(feature_tbl) || nrow(feature_tbl) == 0) {
        stop(
          "El objeto RDS existe, pero no contiene una tabla válida para PANEL_CHOICE = ",
          panel_choice,
          "\nArchivo: ", PANEL_OBJECT_FILE
        )
      }
      
      feature_tbl <- tibble::as_tibble(feature_tbl)
      
      return(list(
        feature_tbl = feature_tbl,
        feature_file_used = PANEL_OBJECT_FILE,
        feature_input_type = "panel_multimodel_rds",
        diagnostic_source = "PANEL_FINAL"
      ))
    }
    
    PANEL_FILE_CSV <- switch(
      panel_choice,
      "FINAL" = file.path(PANEL_TABLES_DIR, "07_PANEL_FINAL_mofa_2of3.csv"),
      "ALL"   = file.path(PANEL_TABLES_DIR, "06_panel_feature_consensus_all.csv"),
      stop("PANEL_CHOICE no reconocido: ", panel_choice)
    )
    
    if (!file.exists(PANEL_FILE_CSV)) {
      stop(
        "No existe la salida del modelamiento multimodelo:\n",
        PANEL_FILE_CSV,
        "\nPrimero corre panel_final_multimodelo.R para ANALYSIS_NAME = ",
        ANALYSIS_NAME,
        "\nMODEL_SOURCE = ", MODEL_SOURCE,
        "\nPANEL_CHOICE = ", PANEL_CHOICE
      )
    }
    
    cat("\nPanel cargado desde CSV del modelamiento multimodelo:\n")
    cat(PANEL_FILE_CSV, "\n")
    
    feature_tbl <- read.csv(PANEL_FILE_CSV, check.names = FALSE) %>%
      tibble::as_tibble()
    
    return(list(
      feature_tbl = feature_tbl,
      feature_file_used = PANEL_FILE_CSV,
      feature_input_type = "panel_multimodel_csv",
      diagnostic_source = "PANEL_FINAL"
    ))
  }
  
  stop("DIAG_TARGET no reconocido: ", diag_target)
}

modelo <- readRDS(MOFA_MODEL_FILE)
rfm <- readRDS(RFM_FILE)

trans_metadata <- rfm$features_metadata

df_train <- readRDS(DF_TRAIN_FILE)
df_test  <- readRDS(DF_TEST_FILE)

diagnostic_loaded <- load_diagnostic_feature_set(
  diag_target = DIAG_TARGET,
  panel_choice = PANEL_CHOICE
)

## Mantengo nombre panel_tbl para no romper el resto del script.
## Pero ahora puede representar:
##   - features finales BRMS/MOFA
##   - panel final multimodelo
panel_tbl <- diagnostic_loaded$feature_tbl
PANEL_FILE <- diagnostic_loaded$feature_file_used
PANEL_INPUT_TYPE <- diagnostic_loaded$feature_input_type
DIAGNOSTIC_SOURCE <- diagnostic_loaded$diagnostic_source

if (!"feature_model" %in% colnames(panel_tbl)) {
  stop(
    "La tabla diagnóstica no tiene columna feature_model: ",
    PANEL_FILE
  )
}

panel_tbl <- panel_tbl %>%
  dplyr::distinct(feature_model, .keep_all = TRUE)

## Objeto compatible con el script antiguo
biomarcadores <- list(
  features = panel_tbl %>%
    dplyr::transmute(
      Feature = feature_model,
      Candidate = "Yes"
    ) %>%
    dplyr::distinct()
)

feats_included <- colnames(df_train)

all_data <- dplyr::bind_rows(df_train, df_test)
## ============================================================
## Añadir factores MOFA solo para diagnóstico
## No se usan como predictores del modelo final
## ============================================================

if (!all(c("Factor1", "Factor2") %in% colnames(all_data))) {
  
  Z_train <- MOFA2::get_factors(
    modelo,
    factors = 1:2,
    as.data.frame = FALSE
  )
  
  if (is.list(Z_train)) {
    Z_train <- Z_train[[1]]
  }
  
  Z_train <- as.matrix(Z_train)
  
  ids_train <- rfm$grupo_train$id
  ids_test  <- rfm$grupo_test$id
  
  ## Asegurar orientación muestras x factores
  if (!all(ids_train %in% rownames(Z_train)) && all(ids_train %in% colnames(Z_train))) {
    Z_train <- t(Z_train)
  }
  
  Z_train <- Z_train[ids_train, 1:2, drop = FALSE]
  colnames(Z_train) <- c("Factor1", "Factor2")
  
  W_proj <- MOFA2::get_expectations(modelo, variable = "W")
  
  if (
    is.list(W_proj) &&
    length(W_proj) == 1 &&
    is.list(W_proj[[1]]) &&
    !is.matrix(W_proj[[1]])
  ) {
    W_proj <- W_proj[[1]]
  }
  
  view_map <- c(
    transcriptomica = "tx",
    proteomica      = "pr",
    metabolomica    = "me",
    clinical        = "cl"
  )
  
  Z_test_list <- lapply(names(view_map), function(v) {
    
    code <- view_map[[v]]
    
    if (!v %in% names(W_proj)) {
      return(NULL)
    }
    
    Wv <- as.matrix(W_proj[[v]])
    Xv <- as.matrix(rfm[[code]]$test)
    
    common_feats <- intersect(rownames(Wv), rownames(Xv))
    
    if (length(common_feats) < 2) {
      return(NULL)
    }
    
    W_sub <- Wv[common_feats, 1:2, drop = FALSE]
    X_sub <- Xv[common_feats, , drop = FALSE]
    
    Zv <- t(MASS::ginv(W_sub) %*% X_sub)
    colnames(Zv) <- c("Factor1", "Factor2")
    rownames(Zv) <- colnames(X_sub)
    
    Zv
  })
  
  Z_test_list <- Z_test_list[!vapply(Z_test_list, is.null, logical(1))]
  
  if (length(Z_test_list) == 0) {
    warning("No se pudo proyectar TEST a factores MOFA. Se usarán NA para Factor1/Factor2 en test.")
    Z_test <- matrix(
      NA_real_,
      nrow = length(ids_test),
      ncol = 2,
      dimnames = list(ids_test, c("Factor1", "Factor2"))
    )
  } else {
    Z_test <- Reduce("+", Z_test_list) / length(Z_test_list)
    Z_test <- Z_test[ids_test, 1:2, drop = FALSE]
  }
  
  ## Insertar factores detrás de y
  df_train <- dplyr::bind_cols(
    df_train[, "y", drop = FALSE],
    as.data.frame(Z_train[rownames(df_train), , drop = FALSE]),
    df_train[, setdiff(colnames(df_train), c("y", "Factor1", "Factor2")), drop = FALSE]
  )
  
  df_test <- dplyr::bind_cols(
    df_test[, "y", drop = FALSE],
    as.data.frame(Z_test[rownames(df_test), , drop = FALSE]),
    df_test[, setdiff(colnames(df_test), c("y", "Factor1", "Factor2")), drop = FALSE]
  )
  
  all_data <- dplyr::bind_rows(df_train, df_test)
}
folder <- DIAG_DIR

dir_create(folder)

cat("\n--- Objetos cargados ---\n")
cat("MOFA:", MOFA_MODEL_FILE, "\n")
cat("RFM:", RFM_FILE, "\n")
cat("df_train:", DF_TRAIN_FILE, "\n")
cat("df_test:", DF_TEST_FILE, "\n")
cat("DIAG_TARGET:", DIAG_TARGET, "\n")
cat("DIAGNOSTIC_SOURCE:", DIAGNOSTIC_SOURCE, "\n")
cat("Feature set usado:", PANEL_FILE, "\n")
cat("Diagnóstico:", DIAG_DIR, "\n")

##===plots_mofa2=====

###===variance explained====

folder_integracion <- "integracion"
dir_create(file.path(folder, folder_integracion))

ps <- plot_variance_explained(modelo, plot_total = T)

df_factors <- ps[[1]]$data %>%
  mutate(
    factor = recode(factor, "Factor1" = "Factor 1", "Factor2" = "Factor 2"),
    view = recode(
      view,
      "transcriptomica" = "Transcriptomics",
      "proteomica"      = "Proteomics",
      "metabolomica"    = "Metabolomics",
      "clinical"        = "Clinical"
    )
  )

df_total <- ps[[2]]$data %>%
  mutate(
    view = recode(
      view,
      "transcriptomica" = "Transcriptomics",
      "proteomica"      = "Proteomics",
      "metabolomica"    = "Metabolomics",
      "clinical"        = "Clinical"
    )
  )

# Definir paleta
omics_colors <- c(
  "Clinical"       = "blue",
  "Transcriptomics" = "red",
  "Proteomics"      = "black",
  "Metabolomics"    = "olivedrab"
)
# p1
p1 <- ggplot(df_factors, aes(x = factor, y = value, fill = view)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = omics_colors) +
  labs(x = "Factor", y = "Explained variance (%)", fill = "Omics view") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(face = "bold", colour = "black"),   # TODO EN NEGRITA Y NEGRO
    axis.text.x = element_text(colour = "black"),           # xticks en negro
    axis.text.y = element_text(colour = "black"),           # yticks en negro
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# p2
p2 <- ggplot(df_total, aes(x = view, y = R2, fill = view)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = omics_colors) +
  labs(x = "Omics view", y = "Total explained variance (%)") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(face = "bold", colour = "black"),   # TODO EN NEGRITA Y NEGRO
    axis.text.x = element_text(colour = "black"),           # xticks en negro
    axis.text.y = element_text(colour = "black"),           # yticks en negro
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

fig <- (p1 + p2) + plot_annotation(tag_levels = "A")

fig

ggsave(
  filename = file.path(folder, folder_integracion, "./Variance_Explained.jpeg"),
  plot = fig,
  width = 13,
  height = 8
)


###====factores=======

Z <- as.data.frame(get_factors(modelo)[[1]])
Z$grupo <- factor(modelo@samples_metadata$grupo)
Z$sample <- modelo@samples_metadata$sample  # etiquetas

p_scatter <- ggplot(Z, aes(
  x = Factor1,
  y = Factor2,
  color = grupo,
  label = sample
)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 3,
                  max.overlaps = 20,
                  show.legend = FALSE) +
  labs(
    title = "MOFA2 Latent Factors",
    x = "Factor 1",
    y = "Factor 2",
    color = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Mostrar
p_scatter

ggsave(
  filename = file.path(folder, folder_integracion, "./MOFA2_latent_factors.jpeg"),
  plot = p_scatter,
  width = 8,
  height = 8
)

###====pesos=======


aux <- plot_top_weights(modelo,
                        view = 1:4,
                        factors = 1:2,
                        nfeatures = Inf)
aux$data$feature

# Definir colores consistentes
omics_colors <- c(
  "Clinical"        = "blue",
  "Transcriptomics" = "red",
  "Proteomics"      = "black",
  "Metabolomics"    = "olivedrab"
)

# --- Preparar df_weights ---
# --- Preparar df_weights ---
sanitize_names <- function(v) {
  v2 <- gsub("[^A-Za-z0-9]+", "_", v)
  v2 <- gsub("__+", "_", v2)
  v2 <- gsub("^_+|_+$", "", v2)
  v2 <- ifelse(grepl("^[0-9]", v2), paste0("X", v2), v2)
  v2 <- make.names(v2, unique = FALSE)
  v2 <- gsub("\\.+", "_", v2)
  v2 <- gsub("^_+|_+$", "", v2)
  v2 <- gsub("__+", "_", v2)
  make.unique(v2, sep = "_")
}

view_code_from_view <- function(view) {
  view <- as.character(view)
  
  dplyr::case_when(
    view == "transcriptomica" ~ "tx",
    view == "proteomica"      ~ "pr",
    view == "metabolomica"    ~ "me",
    view == "clinical"        ~ "cl",
    view == "Transcriptomics" ~ "tx",
    view == "Proteomics"      ~ "pr",
    view == "Metabolomics"    ~ "me",
    view == "Clinical"        ~ "cl",
    TRUE ~ NA_character_
  )
}

view_pretty <- function(view) {
  view <- as.character(view)
  
  dplyr::case_when(
    view == "transcriptomica" ~ "Transcriptomics",
    view == "proteomica"      ~ "Proteomics",
    view == "metabolomica"    ~ "Metabolomics",
    view == "clinical"        ~ "Clinical",
    view == "Transcriptomics" ~ "Transcriptomics",
    view == "Proteomics"      ~ "Proteomics",
    view == "Metabolomics"    ~ "Metabolomics",
    view == "Clinical"        ~ "Clinical",
    TRUE ~ view
  )
}
df_weights <- aux$data %>%
  dplyr::mutate(
    view_original = as.character(view),
    view = view_pretty(view),
    view_code = view_code_from_view(view),
    weight_signed = ifelse(sign == "+", value, -value),
    feature_id = as.character(feature)
  ) %>%
  dplyr::left_join(
    trans_metadata %>%
      dplyr::distinct(EntrezGeneID, .keep_all = TRUE) %>%
      dplyr::transmute(
        EntrezGeneID = as.character(EntrezGeneID),
        GeneSymbol = as.character(GeneSymbol)
      ),
    by = c("feature_id" = "EntrezGeneID")
  ) %>%
  dplyr::mutate(
    label_pretty = dplyr::case_when(
      view == "Transcriptomics" & !is.na(GeneSymbol) & GeneSymbol != "" ~ GeneSymbol,
      TRUE ~ feature_id
    ),
    
    ## Nombre real usado en all_data / df_train / df_test
    feature_model = sanitize_names(paste0(view_code, "_", label_pretty)),
    
    ## Para compatibilidad con el resto del script:
    ## label ahora debe ser el nombre real de columna
    label = feature_model
  )

cat("\n--- Auditoría cruce df_weights vs all_data ---\n")
cat("df_weights antes de filtrar:", dim(df_weights), "\n")
cat("Features en all_data:", length(setdiff(colnames(all_data), c("y", "Factor1", "Factor2"))), "\n")
cat("Cruce feature_model:", length(intersect(df_weights$feature_model, colnames(all_data))), "\n")

df_weights <- df_weights %>%
  dplyr::filter(feature_model %in% colnames(all_data))

cat("df_weights después de filtrar:", dim(df_weights), "\n")

if (nrow(df_weights) == 0) {
  stop("df_weights quedó vacío: no hay cruce entre feature_model y colnames(all_data).")
}
## ============================================================
## PANEL FINAL: vector único para todos los plots diagnósticos
## ============================================================

panel_features <- biomarcadores$features$Feature
panel_features <- intersect(panel_features, colnames(all_data))

cat("\n--- Auditoría panel_features ---\n")
cat("Features en biomarcadores:", nrow(biomarcadores$features), "\n")
cat("Features del panel presentes en all_data:", length(panel_features), "\n")
print(panel_features)

if (length(panel_features) < 2) {
  stop("Muy pocas variables del panel presentes en all_data.")
}

## Datos no integrados SOLO con panel
all_data_panel <- all_data %>%
  dplyr::select(
    dplyr::any_of(c("y", "Factor1", "Factor2")),
    dplyr::all_of(panel_features)
  )

## Pesos MOFA solo del panel
df_weights_panel <- df_weights %>%
  dplyr::filter(feature_model %in% panel_features)

## Anotación única del panel para heatmaps
annotation_col_panel <- df_weights_panel %>%
  dplyr::distinct(feature_model, .keep_all = TRUE) %>%
  dplyr::transmute(
    Label = feature_model,
    Label_pretty = label_pretty,
    View = factor(
      view,
      levels = c("Transcriptomics", "Proteomics", "Metabolomics", "Clinical")
    )
  )

rownames(annotation_col_panel) <- annotation_col_panel$Label

df_weights_unique_panel <- df_weights_panel %>%
  dplyr::distinct(feature_model, view, label_pretty) %>%
  dplyr::rename(
    Feature = feature_model
  )
# ===================== #
# Factor 1
df_factor1 <- df_weights %>%
  filter(factor == "Factor1") %>%
  group_by(feature_model) %>%
  slice_max(order_by = abs(weight_signed), n = 1) %>%
  ungroup() %>%
  arrange(weight_signed) %>%
  mutate(order = row_number())

p1 <- ggplot(df_factor1, aes(x = order, y = weight_signed, fill = view)) +
  geom_col() +
  scale_fill_manual(values = omics_colors) +
  coord_flip() +
  scale_x_continuous(breaks = df_factor1$order, labels = df_factor1$label_pretty) +
  labs(title = "MOFA2 Feature Weights (Factor 1)",
       x = "Feature (ordered)",
       y = "Signed weight",
       fill = "Omics view") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold", size = 12)   # <<< AUMENTADO + NEGRITA
  )

# ===================== #
# Factor 2
df_factor2 <- df_weights %>%
  filter(factor == "Factor2") %>%
  group_by(feature_model) %>%
  slice_max(order_by = abs(weight_signed), n = 1) %>%
  ungroup() %>%
  arrange(weight_signed) %>%
  mutate(order = row_number())

p2 <- ggplot(df_factor2, aes(x = order, y = weight_signed, fill = view)) +
  geom_col() +
  scale_fill_manual(values = omics_colors) +
  coord_flip() +
  scale_x_continuous(breaks = df_factor2$order, labels = df_factor2$label_pretty) +
  labs(title = "MOFA2 Feature Weights (Factor 2)",
       x = "Feature (ordered)",
       y = "Signed weight",
       fill = "Omics view") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold", size = 12)   # <<< MISMO CAMBIO AQUÍ
  )

# ===================== #
# Figura final
fig <- p1 + p2 + plot_annotation(tag_levels = "A")
fig

ggsave(
  filename = file.path(folder, folder_integracion, "./MOFA2_used_weights.jpeg"),
  plot = fig,
  width = 17,
  height = 15
)

## ============================================================
## MOFA weights SOLO para variables del panel final
## ============================================================

df_factor1_panel <- df_weights_panel %>%
  dplyr::filter(factor == "Factor1") %>%
  dplyr::group_by(feature_model) %>%
  dplyr::slice_max(order_by = abs(weight_signed), n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(weight_signed) %>%
  dplyr::mutate(order = dplyr::row_number())

p1_panel <- ggplot(df_factor1_panel, aes(x = order, y = weight_signed, fill = view)) +
  geom_col() +
  scale_fill_manual(values = omics_colors) +
  coord_flip() +
  scale_x_continuous(
    breaks = df_factor1_panel$order,
    labels = df_factor1_panel$label_pretty
  ) +
  labs(
    title = "MOFA2 Feature Weights - Panel biomarkers (Factor 1)",
    x = "Feature",
    y = "Signed weight",
    fill = "Omics view"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold", size = 12)
  )

df_factor2_panel <- df_weights_panel %>%
  dplyr::filter(factor == "Factor2") %>%
  dplyr::group_by(feature_model) %>%
  dplyr::slice_max(order_by = abs(weight_signed), n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(weight_signed) %>%
  dplyr::mutate(order = dplyr::row_number())

p2_panel <- ggplot(df_factor2_panel, aes(x = order, y = weight_signed, fill = view)) +
  geom_col() +
  scale_fill_manual(values = omics_colors) +
  coord_flip() +
  scale_x_continuous(
    breaks = df_factor2_panel$order,
    labels = df_factor2_panel$label_pretty
  ) +
  labs(
    title = "MOFA2 Feature Weights - Panel biomarkers (Factor 2)",
    x = "Feature",
    y = "Signed weight",
    fill = "Omics view"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold", size = 12)
  )

fig_panel_weights <- p1_panel + p2_panel + patchwork::plot_annotation(tag_levels = "A")
fig_panel_weights

ggsave(
  filename = file.path(folder, folder_integracion, "MOFA2_panel_biomarker_weights.jpeg"),
  plot = fig_panel_weights,
  width = 17,
  height = 15
)


##===bivariante=====
# Colores fijos
omics_colors <- c(
  "Clinical"        = "blue",
  "Transcriptomics" = "red",
  "Proteomics"      = "black",
  "Metabolomics"    = "olivedrab"
)

valid_feats <- panel_features

annotation_col <- annotation_col_panel
ann_colors <- list(View = omics_colors)

p1 <- make_corr_heatmap_global(all_data_panel, annotation_col_panel, ann_colors)
p2 <- make_corr_heatmap_group(all_data_panel, "NP", annotation_col_panel, ann_colors)
p3 <- make_corr_heatmap_group(all_data_panel, "SO", annotation_col_panel, ann_colors)
p4 <- make_corr_heatmap_group(all_data_panel, "FD", annotation_col_panel, ann_colors)


folder_correlaciones <- "correlaciones"
dir_create(file.path(folder, folder_correlaciones))

ggsave(
  filename = file.path(
    folder,
    folder_correlaciones,
    "./correlaciones_globales.jpeg"
  ),
  plot = p1,
  width = 10,
  height = 8
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_NP.jpeg"),
  plot = p2,
  width = 10,
  height = 8
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_SO.jpeg"),
  plot = p3,
  width = 10,
  height = 8
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_FD.jpeg"),
  plot = p4,
  width = 10,
  height = 8
)




# --- Colores de ómicas ---
omics_colors <- c(
  "Transcriptomics" = "red",
  "Proteomics"      = "black",
  "Metabolomics"    = "olivedrab",
  "Clinical"        = "blue"
)

# --- Anotación de features ---
valid_feats <- panel_features

annotation_col <- annotation_col_panel
ann_colors <- list(View = omics_colors)

p_global <- make_feature_corr_heatmap(all_data_panel, annotation_col_panel, ann_colors, "Global")
p_NP     <- make_feature_corr_heatmap(all_data_panel, annotation_col_panel, ann_colors, "NP")
p_SO     <- make_feature_corr_heatmap(all_data_panel, annotation_col_panel, ann_colors, "SO")
p_FD     <- make_feature_corr_heatmap(all_data_panel, annotation_col_panel, ann_colors, "FD")


folder_correlaciones <- "correlaciones_bivariante"
dir_create(file.path(folder, folder_correlaciones))

ggsave(
  filename = file.path(
    folder,
    folder_correlaciones,
    "./correlaciones_globales_sin_integrar.jpeg"
  ),
  plot = p_global,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_NP_sin_integrar.jpeg"),
  plot = p_NP,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_SO_sin_integrar.jpeg"),
  plot = p_SO,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_FD_sin_integrar.jpeg"),
  plot = p_FD,
  width = 10,
  height = 10
)


## ============================================================
## Datos TEST para proyección diagnóstica al espacio MOFA
## ============================================================

tst_Data <- rfm[1:4]
tst_Data <- lapply(tst_Data, function(x) x$test)
tst_Data <- t(Reduce("rbind", tst_Data))

Z <- MOFA2::get_expectations(modelo, variable = "Z")
W <- MOFA2::get_expectations(modelo, variable = "W")

# --- Función para proyectar nuevos datos al espacio latente de MOFA2 ---
project_to_mofa <- function(tst_Data, W) {
  # tst_Data: matriz (features × muestras)
  # W: lista de pesos de MOFA2, uno por vista
  
  # asegurar que las features estén en filas
  if (nrow(tst_Data) < ncol(tst_Data)) {
    tst_Data <- t(tst_Data)
  }
  
  # Separar tst_Data en lista de vistas usando rownames(W)
  tst_list <- lapply(W, function(Wv) {
    feats <- intersect(rownames(Wv), rownames(tst_Data))
    tst_Data[feats, , drop = FALSE]
  })
  
  # Proyección vista por vista
  Z_estimates <- lapply(names(W), function(vista) {
    X <- tst_list[[vista]]
    Wv <- W[[vista]]
    
    feats <- intersect(rownames(Wv), rownames(X))
    X <- X[feats, , drop = FALSE]
    Wv <- Wv[feats, , drop = FALSE]
    
    # pseudoinversa de Wv
    pinv <- solve(t(Wv) %*% Wv) %*% t(Wv)
    Zv <- t(pinv %*% X)   # muestras × factores
    Zv
  })
  
  # Combinar proyecciones promediando entre vistas
  Z_combined <- Reduce("+", Z_estimates) / length(Z_estimates)
  return(Z_combined)
}

# --- Ejemplo de uso ---
# tst_Data ya cargado (10 × 504)
# W ya calculado con get_expectations(modelo, "W")

W_ <- Reduce("rbind",W)

Z_new <- project_to_mofa(tst_Data, W)

dim(Z_new)  # debería ser 10 muestras × 2 factores
head(Z_new)

R_new <- Z_new %*% t(W_)


R <- MOFA2::predict(modelo)
R <- lapply(R, function(x) x[[1]])
R <- t(Reduce("rbind",R))
R <- as.data.frame(R)


all_R <- as.data.frame(rbind(R,R_new))
# Crear diccionario EntrezID -> GeneSymbol (solo transcriptómica)
map_entrez2symbol <- trans_metadata %>%
  filter(!is.na(EntrezGeneID), !is.na(GeneSymbol)) %>%
  distinct(EntrezGeneID, GeneSymbol)

# Aseguramos que EntrezGeneID sea texto, igual que tus colnames
map_entrez2symbol <- map_entrez2symbol %>%
  mutate(EntrezGeneID = as.character(EntrezGeneID))

# Reemplazar los nombres de las columnas en R
new_names <- colnames(all_R)

# Solo cambiamos donde hay match
new_names <- ifelse(new_names %in% map_entrez2symbol$EntrezGeneID,
                    map_entrez2symbol$GeneSymbol[match(new_names, map_entrez2symbol$EntrezGeneID)],
                    new_names)

colnames(all_R) <- new_names

## ============================================================
## Mapeo seguro: features del modelo -> columnas reconstruidas MOFA
## No usar fuzzy matching aquí
## ============================================================

train_feats <- setdiff(colnames(df_train), c("y", "Factor1", "Factor2"))

map_feats <- df_weights %>%
  dplyr::distinct(feature_model, label_pretty, feature_id, view) %>%
  dplyr::filter(feature_model %in% train_feats) %>%
  dplyr::mutate(
    R_col_candidate_1 = label_pretty,
    R_col_candidate_2 = feature_id,
    R_col = dplyr::case_when(
      R_col_candidate_1 %in% colnames(all_R) ~ R_col_candidate_1,
      R_col_candidate_2 %in% colnames(all_R) ~ R_col_candidate_2,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::arrange(match(feature_model, train_feats))

cat("\n--- Auditoría mapeo all_R ---\n")
cat("Features train:", length(train_feats), "\n")
cat("Features mapeadas en df_weights:", nrow(map_feats), "\n")
cat("Features con R_col válido:", sum(!is.na(map_feats$R_col)), "\n")

cat("\n--- Features sin match en all_R ---\n")
print(map_feats %>% dplyr::filter(is.na(R_col)))

if (any(is.na(map_feats$R_col))) {
  warning("Hay features sin match en all_R. Se excluirán solo del análisis integrado reconstruido.")
}

map_feats_ok <- map_feats %>%
  dplyr::filter(!is.na(R_col))

R_sub <- all_R[, map_feats_ok$R_col, drop = FALSE]

## Renombrar al formato del modelo: tx_..., pr_..., me_..., cl_...
colnames(R_sub) <- map_feats_ok$feature_model

cat("\n--- Dim R_sub reconstruido ---\n")
print(dim(R_sub))

cat("\n--- Primeras columnas R_sub ---\n")
print(head(colnames(R_sub), 30))
### biomarcadores 65 !!!!!

biomarkers <- biomarcadores
## ============================================================
## Subset al panel final seleccionado
## Ahora R_sub ya tiene nombres feature_model
## ============================================================

panel_features_integrado <- intersect(
  biomarkers$features$Feature,
  colnames(R_sub)
)

cat("\n--- Panel vs R_sub reconstruido ---\n")
cat("Features panel:", nrow(biomarkers$features), "\n")
cat("Features panel presentes en R_sub:", length(panel_features_integrado), "\n")

cat("\n--- Features del panel NO presentes en R_sub ---\n")
print(setdiff(biomarkers$features$Feature, colnames(R_sub)))

if (length(panel_features_integrado) < 2) {
  stop("Muy pocas features del panel presentes en R_sub reconstruido.")
}

R_sub <- R_sub[, panel_features_integrado, drop = FALSE]

cat("\n--- Dim R_sub tras filtrar por panel ---\n")
print(dim(R_sub))

y_all <- sub("_.*", "", rownames(R_sub))
R_sub$y <-y_all
# --- Ejecutar ---
p_global <- make_feature_corr_heatmap(R_sub, annotation_col_panel, ann_colors, "Global")
p_NP     <- make_feature_corr_heatmap(R_sub, annotation_col_panel, ann_colors, "NP")
p_SO     <- make_feature_corr_heatmap(R_sub, annotation_col_panel, ann_colors, "SO")
p_FD     <- make_feature_corr_heatmap(R_sub, annotation_col_panel, ann_colors, "FD") ## must have more than >4 observation its has 3

cor_NP <- psych::corr.test(R_sub[R_sub$y=="NP",-ncol(R_sub)],method = "spearman",adjust = "BH")
write.csv(cor_NP$r, out_table("cor_NP_integrado_r.csv"))
write.csv(cor_NP$p, out_table("cor_NP_integrado_p.csv"))

cor_SO <- psych::corr.test(R_sub[R_sub$y=="SO",-ncol(R_sub)],method = "spearman",adjust = "BH")
write.csv(cor_SO$r, out_table("cor_SO_integrado_r.csv"))
write.csv(cor_SO$p, out_table("cor_SO_integrado_p.csv"))

cor_FD <- psych::corr.test(R_sub[R_sub$y=="FD",-ncol(R_sub)],method = "spearman",adjust = "BH")
write.csv(cor_FD$r, out_table("cor_FD_integrado_r.csv"))
write.csv(cor_FD$p, out_table("cor_FD_integrado_p.csv"))


folder_correlaciones <- "correlaciones_bivariante_integrado"
dir_create(file.path(folder, folder_correlaciones))

ggsave(
  filename = file.path(
    folder,
    folder_correlaciones,
    "./correlaciones_globales_integrado.jpeg"
  ),
  plot = p_global,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_NP_integrado.jpeg"),
  plot = p_NP,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_SO_integrado.jpeg"),
  plot = p_SO,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_FD_integrado.jpeg"),
  plot = p_FD,
  width = 10,
  height = 10
)



##====multivariante ======

## ============================================================
## PCA DIAGNÓSTICO
## Ambos PCA usan exactamente las variables del panel:
## 1) Sin integrar  : all_data_panel
## 2) Integrado     : R_sub reconstruido por MOFA
## ============================================================

folder_multivariante <- "multivariante"
dir_create(file.path(folder, folder_multivariante))

df_weights_unique <- df_weights_unique_panel

panel_features_pca <- intersect(panel_features, colnames(all_data_panel))
panel_features_pca <- intersect(panel_features_pca, colnames(R_sub))

cat("\n--- Auditoría PCA panel ---\n")
cat("Features panel:", length(panel_features), "\n")
cat("Features panel para PCA sin/integrado:", length(panel_features_pca), "\n")
print(panel_features_pca)

if (length(panel_features_pca) < 2) {
  stop("Muy pocas variables del panel para PCA.")
}

## ------------------------------------------------------------
## 1) PCA SIN INTEGRAR: valores originales del panel
## ------------------------------------------------------------

X_pca_sin_integrar <- all_data_panel[, panel_features_pca, drop = FALSE]

pcx_sin_integrar <- prcomp(
  X_pca_sin_integrar,
  center = TRUE,
  scale. = TRUE
)

p_loading_sin_integrar <- plot_pca_feats(
  pcx_sin_integrar,
  df_weights_unique,
  panel_features_pca,
  axes = c(1, 2)
) +
  ggtitle("PCA loadings - Panel biomarkers, non-integrated data")

varianza_sin <- round(
  100 * (pcx_sin_integrar$sdev^2) / sum(pcx_sin_integrar$sdev^2),
  2
)

scores_sin <- data.frame(
  pcx_sin_integrar$x[, 1:2, drop = FALSE],
  Group = all_data_panel$y,
  name = rownames(all_data_panel)
)

p_score_sin_integrar <- ggplot(
  scores_sin,
  aes(PC1, PC2, color = Group, label = name)
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(max.overlaps = 50, size = 3) +
  labs(
    title = "PCA score plot - Panel biomarkers, non-integrated data",
    x = paste0("PC1 (", varianza_sin[1], "%)"),
    y = paste0("PC2 (", varianza_sin[2], "%)"),
    color = "Group"
  ) +
  scale_color_manual(values = c("NP" = "black", "SO" = "blue", "FD" = "red")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  )

## ------------------------------------------------------------
## 2) PCA INTEGRADO: reconstrucción MOFA del panel
## ------------------------------------------------------------

X_pca_integrado <- R_sub[, panel_features_pca, drop = FALSE]

pcx_integrado <- prcomp(
  X_pca_integrado,
  center = TRUE,
  scale. = TRUE
)

p_loading_integrado <- plot_pca_feats(
  pcx_integrado,
  df_weights_unique,
  panel_features_pca,
  axes = c(1, 2)
) +
  ggtitle("PCA loadings - Panel biomarkers, MOFA-reconstructed data")

varianza_int <- round(
  100 * (pcx_integrado$sdev^2) / sum(pcx_integrado$sdev^2),
  2
)

scores_int <- data.frame(
  pcx_integrado$x[, 1:2, drop = FALSE],
  Group = R_sub$y,
  name = rownames(R_sub)
)

p_score_integrado <- ggplot(
  scores_int,
  aes(PC1, PC2, color = Group, label = name)
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(max.overlaps = 50, size = 3) +
  labs(
    title = "PCA score plot - Panel biomarkers, MOFA-reconstructed data",
    x = paste0("PC1 (", varianza_int[1], "%)"),
    y = paste0("PC2 (", varianza_int[2], "%)"),
    color = "Group"
  ) +
  scale_color_manual(values = c("NP" = "black", "SO" = "blue", "FD" = "red")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  )

## ------------------------------------------------------------
## Guardar PCA
## ------------------------------------------------------------

ggsave(
  filename = file.path(folder, folder_multivariante, "pca_panel_loadings_sin_integrar.jpeg"),
  plot = p_loading_sin_integrar,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_multivariante, "pca_panel_loadings_integrado.jpeg"),
  plot = p_loading_integrado,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_multivariante, "pca_panel_score_sin_integrar.jpeg"),
  plot = p_score_sin_integrar,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_multivariante, "pca_panel_score_integrado.jpeg"),
  plot = p_score_integrado,
  width = 10,
  height = 10
)

p_loading_sin_integrar
p_loading_integrado
p_score_sin_integrar
p_score_integrado
##===univariante=====

###====integrado=====


R_sub$y <- as.factor(R_sub$y)

FD_NP <- lapply(colnames(R_sub)[-ncol(R_sub)], function(col) {
  csa <- contrast_t_test(
    R_sub,
    factor_ = "y",
    response = col,
    contrast_expr = "FD-NP"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})

SO_NP <- lapply(colnames(R_sub)[-ncol(R_sub)], function(col) {
  csa <- contrast_t_test(
    R_sub,
    factor_ = "y",
    response = col,
    contrast_expr = "SO-NP"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})


FD_SO <- lapply(colnames(R_sub)[-ncol(R_sub)], function(col) {
  csa <- contrast_t_test(
    R_sub,
    factor_ = "y",
    response = col,
    contrast_expr = "SO-FD"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})

FD_NP <- Reduce("rbind",FD_NP)
rownames(FD_NP) <- colnames(R_sub)[-ncol(R_sub)]

SO_NP <- Reduce("rbind",SO_NP)
rownames(SO_NP) <- colnames(R_sub)[-ncol(R_sub)]

FD_SO <- Reduce("rbind",FD_SO)
rownames(FD_SO) <- colnames(R_sub)[-ncol(R_sub)]


statistic <-  data.frame(FD_NP = FD_NP[,1],
                         SO_NP = SO_NP[,1],
                         FD_SO = FD_SO[,1])

pval <-  data.frame(FD_NP = FD_NP[,2],
                         SO_NP = SO_NP[,2],
                         FD_SO = FD_SO[,2])

p_adj <- apply(pval, 2, function(col) {
  out <- rep(NA_real_, length(col))
  ok <- is.finite(col)
  out[ok] <- p.adjust(col[ok], method = "BH")
  out
})

rownames(p_adj) <- rownames(pval)


p1 <- plot_heatmap_stats_ordered(
  statistic = statistic,
  p_adj = p_adj,
  df_weights = df_weights_unique_panel,
  alpha = 0.05
)
folder_correlaciones <- "univariante_integrado"

dir_create(file.path(folder, folder_correlaciones))

ggsave(
  filename = file.path(folder, folder_correlaciones, "./univariante_integrado.jpeg"),
  plot = p1,
  width = 15,
  height = 10
)


###====integrado=====

all_data_ <- all_data_panel %>%
  dplyr::select(-dplyr::any_of(c("Factor1", "Factor2")))


FD_NP <- lapply(colnames(all_data_)[-1], function(col) {
  csa <- contrast_t_test(
    all_data_,
    factor_ = "y",
    response = col,
    contrast_expr = "FD-NP"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})

SO_NP <- lapply(colnames(all_data_)[-1], function(col) {
  csa <- contrast_t_test(
    all_data_,
    factor_ = "y",
    response = col,
    contrast_expr = "SO-NP"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})


FD_SO <- lapply(colnames(all_data_)[-1], function(col) {
  csa <- contrast_t_test(
    all_data_,
    factor_ = "y",
    response = col,
    contrast_expr = "FD-SO"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})

FD_NP <- Reduce("rbind",FD_NP)
rownames(FD_NP) <- colnames(all_data_)[-1]

SO_NP <- Reduce("rbind",SO_NP)
rownames(SO_NP) <- colnames(all_data_)[-1]

FD_SO <- Reduce("rbind",FD_SO)
rownames(FD_SO) <- colnames(all_data_)[-1]


statistic <-  data.frame(FD_NP = FD_NP[,1],
                         SO_NP = SO_NP[,1],
                         FD_SO = FD_SO[,1])

pval <-  data.frame(FD_NP = FD_NP[,2],
                    SO_NP = SO_NP[,2],
                    FD_SO = FD_SO[,2])

p_adj <- apply(pval, 2, function(col) {
  out <- rep(NA_real_, length(col))
  ok <- is.finite(col)
  out[ok] <- p.adjust(col[ok], method = "BH")
  out
})

rownames(p_adj) <- rownames(pval)




p1 <- plot_heatmap_stats_ordered(statistic = statistic,
                           p_adj = p_adj ,
                           df_weights = df_weights_unique_panel,
                           alpha = 0.05)


ggsave(
  filename = file.path(folder, folder_correlaciones, "./univariante_sin_integrado.jpeg"),
  plot = p1,
  width = 15,
  height = 10
)
cat("\n============================================================\n")
cat("DIAGNÓSTICO TERMINADO\n")
cat("Resultados guardados en:\n")
cat("Diagnóstico:", DIAG_DIR, "\n")
cat("Tablas     :", DIAG_TABLES_DIR, "\n")
cat("RDS        :", DIAG_RDS_DIR, "\n")
cat("Panel usado:", PANEL_FILE, "\n")
cat("============================================================\n")
## ============================================================
## GUARDAR OBJETOS DIAGNÓSTICOS
## ============================================================

diagnostic_config <- list(
  analysis_name = ANALYSIS_NAME,
  outdir_base = OUTDIR_BASE,
  diag_target = DIAG_TARGET,
  diagnostic_source = DIAGNOSTIC_SOURCE,
  panel_choice = PANEL_CHOICE,
  model_source = MODEL_SOURCE,
  model_dir_name = MODEL_DIR_NAME,
  panel_dir_name = PANEL_DIR_NAME,
  analysis_dir = ANALYSIS_DIR,
  diag_dir = DIAG_DIR,
  
  rfm_file = RFM_FILE,
  mofa_model_file = MOFA_MODEL_FILE,
  
  modelamiento_rds_dir = MODELAMIENTO_RDS_DIR,
  modelamiento_tables_dir = MODELAMIENTO_TABLES_DIR,
  df_train_file = DF_TRAIN_FILE,
  df_test_file = DF_TEST_FILE,
  
  bayes_features_rds_file = BAYES_FEATURES_RDS_FILE,
  bayes_features_csv_file = BAYES_FEATURES_CSV_FILE,
  bayes_features_object_file = BAYES_FEATURES_OBJECT_FILE,
  
  panel_tables_dir = PANEL_TABLES_DIR,
  panel_rds_dir = PANEL_RDS_DIR,
  panel_object_file = PANEL_OBJECT_FILE,
  
  feature_file_used = PANEL_FILE,
  feature_input_type = PANEL_INPUT_TYPE
)

saveRDS(
  diagnostic_config,
  out_rds("diagnostic_config.rds")
)

saveRDS(
  list(
    diagnostic_config = diagnostic_config,
    panel_tbl = panel_tbl,
    biomarcadores = biomarcadores,
    df_train = df_train,
    df_test = df_test,
    all_data = all_data,
    df_weights = df_weights,
    df_weights_unique = df_weights_unique
  ),
  out_rds("diagnostic_objects.rds")
)

write.csv(
  panel_tbl,
  out_table(paste0("features_usadas_para_diagnostico_", tolower(DIAG_TARGET), ".csv")),
  row.names = FALSE
)

## ============================================================
## Correlación entre factores latentes MOFA
## plot_factor_cor() no devuelve ggplot; devuelve matriz
## Por eso se guarda abriendo dispositivo gráfico
## ============================================================

factor_cor_file <- file.path(
  folder,
  folder_integracion,
  "MOFA2_factor_correlation.jpeg"
)

grDevices::jpeg(
  filename = factor_cor_file,
  width = 8,
  height = 7,
  units = "in",
  res = 300
)

factor_cor_mat <- MOFA2::plot_factor_cor(modelo)

grDevices::dev.off()

## Guardar también la matriz numérica
write.csv(
  factor_cor_mat,
  out_table("MOFA2_factor_correlation_matrix.csv"),
  row.names = TRUE
)

cat("\nFactor correlation plot guardado en:\n")
cat(factor_cor_file, "\n")
cat("Matriz guardada en:\n")
cat(out_table("MOFA2_factor_correlation_matrix.csv"), "\n")
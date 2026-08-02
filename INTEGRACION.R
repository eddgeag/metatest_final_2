library(MOFA2)

# ============================================================
# CONFIGURACIÓN PARAMETRIZABLE DEL ANÁLISIS
# Debe usar el mismo ANALYSIS_NAME que preprocess_data.R
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

ANALYSIS_NAME <- get_env_chr("ANALYSIS_NAME", required = TRUE)
OUTDIR_BASE   <- get_env_chr("OUTDIR_BASE", ".")

OUTPUT_DIR <- file.path(OUTDIR_BASE, ANALYSIS_NAME)
RFM_FILE   <- file.path(
  OUTPUT_DIR,
  paste0("ready_for_modeling_", ANALYSIS_NAME, ".rds")
)
INTEGRATION_DIR <- file.path(OUTPUT_DIR, "integracion")
MOFA_HDF5_DIR   <- file.path(INTEGRATION_DIR, "modelos_hdf5")
MOFA_RDS_DIR    <- file.path(INTEGRATION_DIR, "rds")
MOFA_PLOTS_DIR  <- file.path(INTEGRATION_DIR, "plots")

dir.create(OUTPUT_DIR,      recursive = TRUE, showWarnings = FALSE)
dir.create(INTEGRATION_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MOFA_HDF5_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(MOFA_RDS_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(MOFA_PLOTS_DIR,  recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("ANÁLISIS MOFA:", ANALYSIS_NAME, "\n")
cat("OUTDIR_BASE:", OUTDIR_BASE, "\n")
cat("OUTPUT_DIR:", OUTPUT_DIR, "\n")
cat("RFM_FILE:", RFM_FILE, "\n")
cat("INTEGRATION_DIR:", INTEGRATION_DIR, "\n")
cat("MOFA_HDF5_DIR:", MOFA_HDF5_DIR, "\n")
cat("MOFA_RDS_DIR:", MOFA_RDS_DIR, "\n")
cat("MOFA_PLOTS_DIR:", MOFA_PLOTS_DIR, "\n")
cat("============================================================\n")
## ------------------------------------------------------------
## 0) Helpers
## ------------------------------------------------------------

.get_train_ids <- function(rfm) {
  for (nm in c("cl", "tx", "pr", "me")) {
    if (!is.null(rfm[[nm]]) && !is.null(rfm[[nm]]$train)) {
      ids <- colnames(rfm[[nm]]$train)
      if (is.null(ids)) stop(paste0("Faltan colnames en ", nm, "$train"))
      return(ids)
    }
  }
  stop("No hay vistas con matriz $train y colnames definidos.")
}

.stack_view_train <- function(view) {
  Xtr <- view$train
  if (is.null(colnames(Xtr))) stop("La matriz train no tiene colnames, es decir sujetos.")
  return(Xtr)
}

.build_metadata <- function(rfm, extra_cols = NULL) {
  sample_ids <- .get_train_ids(rfm)
  
  if (!is.null(rfm$grupo_df) && !is.null(rfm$grupo_df$train)) {
    gdf <- as.data.frame(rfm$grupo_df$train, stringsAsFactors = FALSE)
    stopifnot(all(c("id", "grupo") %in% colnames(gdf)))
    
    rownames(gdf) <- gdf$id
    gdf <- gdf[sample_ids, , drop = FALSE]
    
    fac <- factor(gdf$grupo)
    names(fac) <- gdf$id
  } else {
    fac <- factor(sub("^([^_]+)_.*$", "\\1", sample_ids))
    names(fac) <- sample_ids
  }
  
  meta <- data.frame(
    sample = sample_ids,
    grupo  = fac,
    stringsAsFactors = FALSE,
    row.names = sample_ids
  )
  
  if (!is.null(extra_cols)) {
    for (nm in names(extra_cols)) {
      v <- extra_cols[[nm]]
      
      if (!is.null(names(v))) {
        v <- v[sample_ids]
      } else if (length(v) == length(sample_ids)) {
        names(v) <- sample_ids
      } else {
        stop(paste0("extra_metadata '", nm, "' no alinea con muestras."))
      }
      
      meta[[nm]] <- v
    }
  }
  
  colnames(meta) <- gsub("[^[:alnum:]]", "_", colnames(meta))
  colnames(meta) <- gsub("_+", "_", colnames(meta))
  colnames(meta) <- sub("_$", "", colnames(meta))
  
  return(meta)
}

.save_mofa_plot <- function(plot_expr, filename, plot_dir,
                            width = 8, height = 6, dpi = 300) {
  
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  outfile <- file.path(plot_dir, filename)
  
  png(
    filename = outfile,
    width = width,
    height = height,
    units = "in",
    res = dpi
  )
  
  p <- plot_expr
  
  if (!is.null(p)) {
    print(p)
  }
  
  dev.off()
  
  invisible(outfile)
}



## ------------------------------------------------------------
## 1) Función central: ahora guarda por factor + semilla
## ------------------------------------------------------------

mofa_componentes_original <- function(ncomp,
                                      semilla,
                                      directorio_modelo,
                                      mofa.obj,
                                      sample_metadata_) {
  
  samples_metadata(mofa.obj) <- sample_metadata_
  
  if (!dir.exists(directorio_modelo)) {
    dir.create(directorio_modelo, recursive = TRUE)
  }
  
  data_opts <- get_default_data_options(mofa.obj)
  data_opts$scale_views <- FALSE
  
  model_opts <- get_default_model_options(mofa.obj)
  model_opts$num_factors <- ncomp
  
  train_opts <- get_default_training_options(mofa.obj)
  train_opts$seed <- semilla
  train_opts$convergence_mode <- "slow"
  
  MOFAobject <- prepare_mofa(
    object = mofa.obj,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  outfile <- file.path(
    directorio_modelo,
    paste0("modelo_f", ncomp, "_seed", semilla, ".hdf5")
  )
  
  if (file.exists(outfile)) {
    message("Eliminando archivo previo: ", outfile)
    file.remove(outfile)
  }
  
  gc()
  reticulate::py_run_string("import gc; gc.collect()")
  
  run_mofa(
    object = MOFAobject,
    outfile = outfile,
    use_basilisk = TRUE,
    save_data = FALSE
  )
}

## ------------------------------------------------------------
## 2) Construcción de vistas SOLO con train
## ------------------------------------------------------------

build_mofa_from_rfm <- function(read_for_modeling) {
  
  vistas_raw <- list(
    transcriptomica = read_for_modeling$tx,
    proteomica      = read_for_modeling$pr,
    metabolomica    = read_for_modeling$me,
    clinical        = read_for_modeling$cl
  )
  
  vistas_raw <- vistas_raw[!vapply(vistas_raw, is.null, logical(1))]
  
  omicas <- lapply(vistas_raw, .stack_view_train)
  
  omicas_MOFA <- create_mofa_from_matrix(omicas)
  
  list(
    omicas = omicas,
    mofa = omicas_MOFA
  )
}

## ------------------------------------------------------------
## 3) Pipeline: ahora factores x semillas
## ------------------------------------------------------------

run_mofa_grid_from_rfm <- function(read_for_modeling,
                                   factores = 2:6,
                                   semillas = c(23355),
                                   directorio_modelo = MOFA_HDF5_DIR,
                                   extra_metadata = NULL,
                                   seleccionar = TRUE,
                                   plot = TRUE,
                                   save_dir = MOFA_RDS_DIR,
                                   plot_dir = MOFA_PLOTS_DIR,
                                   analysis_name = ANALYSIS_NAME) {
  
  built <- build_mofa_from_rfm(read_for_modeling)
  omicas_MOFA <- built$mofa
  
  aux <- .build_metadata(read_for_modeling, extra_cols = extra_metadata)
  
  common_ids <- colnames(built$omicas[[1]])
  
  if (!identical(rownames(aux), common_ids)) {
    stop("rownames(metadata) deben coincidir con colnames de las ómicas train.")
  }
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  saveRDS(built$omicas, file.path(save_dir, "omicas_train.rds"))
  saveRDS(aux, file.path(save_dir, "metadata_train.rds"))
  
  run_config_mofa <- data.frame(
    analysis_name = analysis_name,
    outdir_base = OUTDIR_BASE,
    output_dir = OUTPUT_DIR,
    input_rfm = RFM_FILE,
    integration_dir = INTEGRATION_DIR,
    directorio_modelo = directorio_modelo,
    save_dir = save_dir,
    plot_dir = plot_dir,
    factores = paste(factores, collapse = ";"),
    semillas = paste(semillas, collapse = ";"),
    seleccionar = seleccionar,
    plot = plot,
    stringsAsFactors = FALSE
  )
  write.csv(
    run_config_mofa,
    file.path(save_dir, "run_config_mofa.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    run_config_mofa,
    file.path(save_dir, "run_config_mofa.rds")
  )
  
  grid <- expand.grid(
    ncomp = factores,
    semilla = semillas,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  lista_modelos <- vector("list", nrow(grid))
  
  for (i in seq_len(nrow(grid))) {
    
    ncomp_i <- grid$ncomp[i]
    seed_i  <- grid$semilla[i]
    
    message("Entrenando MOFA | factores = ", ncomp_i, " | seed = ", seed_i)
    
    lista_modelos[[i]] <- mofa_componentes_original(
      ncomp = ncomp_i,
      semilla = seed_i,
      directorio_modelo = directorio_modelo,
      mofa.obj = omicas_MOFA,
      sample_metadata_ = aux
    )
    
    names(lista_modelos)[i] <- paste0("f", ncomp_i, "_seed", seed_i)
  }
  
  names(lista_modelos) <- paste0(
    "f", grid$ncomp,
    "_seed", grid$semilla
  )
  
  MOFA2::compare_elbo(lista_modelos)
  
  modelo <- if (seleccionar) {
    select_model(lista_modelos, plot = plot)
  } else {
    lista_modelos[[length(lista_modelos)]]
  }
  
  if (plot) {
    
    .save_mofa_plot(
      plot_factor_cor(modelo),
      filename = "MOFA_factor_cor.png",
      plot_dir = plot_dir
    )
    
    .save_mofa_plot(
      plot_variance_explained(modelo),
      filename = "MOFA_variance_explained.png",
      plot_dir = plot_dir
    )
    
    .save_mofa_plot(
      plot_factors(modelo, color_by = "grupo"),
      filename = "MOFA_latent_factors_by_grupo.png",
      plot_dir = plot_dir
    )
  }
  saveRDS(
    lista_modelos,
    file.path(save_dir, paste0("lista_modelos_mofa_grid_", analysis_name, ".rds"))
  )
  
  saveRDS(
    modelo,
    file.path(save_dir, paste0("modelo_mofa_seleccionado_", analysis_name, ".rds"))
  )  
  invisible(list(
    modelos = lista_modelos,
    modelo_seleccionado = modelo,
    grid = grid,
    metadata = samples_metadata(modelo)
  ))
}

# ------------------------------------------------------------
# 4) Uso
# ------------------------------------------------------------
# ------------------------------------------------------------
# 4) Uso parametrizado
# ------------------------------------------------------------

if (!file.exists(RFM_FILE)) {
  stop(
    "No existe el archivo de entrada: ", RFM_FILE,
    "\nPrimero corre el script de preprocesamiento para ANALYSIS_NAME = ", ANALYSIS_NAME
  )
}

rfm <- readRDS(RFM_FILE)

out <- run_mofa_grid_from_rfm(
  read_for_modeling = rfm,
  factores = 2:4,
  semillas = c(23355, 101, 202, 303, 404, 2023, 400, 500, 10000, 2000),
  directorio_modelo = MOFA_HDF5_DIR,
  extra_metadata = NULL,
  seleccionar = TRUE,
  plot = TRUE,
  save_dir = MOFA_RDS_DIR,
  plot_dir = MOFA_PLOTS_DIR,
  analysis_name = ANALYSIS_NAME
)


# ------------------------------------------------------------
# 5) Recarga parametrizada de modelos HDF5 y selección final
# ------------------------------------------------------------

modelos_hdf5 <- list.files(
  MOFA_HDF5_DIR,
  full.names = TRUE,
  pattern = "\\.hdf5$"
)

if (!length(modelos_hdf5)) {
  stop("No se encontraron modelos .hdf5 en: ", MOFA_HDF5_DIR)
}

modelos <- lapply(modelos_hdf5, load_model)

names(modelos) <- sub("\\.hdf5$", "", basename(modelos_hdf5))

.save_mofa_plot(
  MOFA2::compare_elbo(modelos),
  filename = "MOFA_compare_elbo.png",
  plot_dir = MOFA_PLOTS_DIR
)

modelo_final <- select_model(modelos, plot = TRUE)

saveRDS(
  modelos,
  file.path(MOFA_RDS_DIR, paste0("modelos_hdf5_cargados_", ANALYSIS_NAME, ".rds"))
)

saveRDS(
  modelo_final,
  file.path(MOFA_RDS_DIR, paste0("modelo_mofa_final_", ANALYSIS_NAME, ".rds"))
)


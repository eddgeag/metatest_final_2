#!/usr/bin/env Rscript
# METATEST2 — script 30b v1.4.1
# Compara curvas ROC one-vs-rest de tres paneles de un MISMO escenario:
# CORE18_47, FINAL28 y MOFA_OTSU. No ajusta modelos ni recalcula probabilidades.
# Se puede ejecutar independientemente, despues del script 30 v1.4.0 o v1.4.1.

options(stringsAsFactors = FALSE)
r30b_families <- c("P1", "P2", "P3")
r30b_classes <- c("NP", "FD", "SO")
r30b_panels <- c("CORE18_47", "FINAL28", "MOFA_OTSU")
r30b_colors <- c(CORE18_47 = "#9171BD", FINAL28 = "#4E79A7", MOFA_OTSU = "#E15759")
r30b_labels <- c(CORE18_47 = "Consenso 47 (script 18)",
                FINAL28 = "Panel final (script 28)",
                MOFA_OTSU = "Pesos MOFA2 (Otsu)")
r30b_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

r30b_choose_source <- function() {
  explicit <- Sys.getenv("MT2_30_ROC_SOURCE_DIR", "")
  if (nzchar(explicit)) {
    r30b_assert(dir.exists(path.expand(explicit)), paste("Carpeta no encontrada:", explicit))
    return(normalizePath(path.expand(explicit), mustWork = TRUE))
  }
  # Si se invoca desde el script 30, reutiliza su carpeta actual.
  if (exists("OUT", inherits = TRUE)) {
    current <- get("OUT", inherits = TRUE)
    if (is.character(current) && length(current) == 1L && dir.exists(current))
      return(normalizePath(current, mustWork = TRUE))
  }
  args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  r30b_assert(length(args) > 0L,
    "Indica MT2_30_ROC_SOURCE_DIR con la carpeta del script 30 ya ejecutado")
  script <- normalizePath(sub("^--file=", "", args[[1L]]), mustWork = TRUE)
  repo <- dirname(dirname(dirname(script)))
  run <- path.expand(Sys.getenv("METATEST2_RUN_ROOT", file.path(
    dirname(repo), "results_cedia", "metatest2_integral_66_20261027T013422Z")))
  sel <- path.expand(Sys.getenv("METATEST2_SELECTION_RUN", file.path(
    run, "comparacion_modelamientos", "upstream_comun", "downstream",
    "seleccion_mejor_escenario_train", "20261029T210729Z_pid2405301")))
  root <- file.path(sel, "SCRIPT30_LOGISTICA_MOFA")
  r30b_assert(dir.exists(root), paste("No existe:", root))
  dirs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[file.exists(file.path(dirs, "tables", "23_ROC_CURVE_POINTS_OVR.csv")) &
                 file.exists(file.path(dirs, "tables", "25_ROC_PR_OVR_BY_CLASS.csv")) &
                 file.exists(file.path(dirs, "tables", "10_MODEL_STATUS.csv"))]
  dirs <- dirs[order(file.info(dirs)$mtime, decreasing = TRUE)]
  for (d in dirs) {
    st <- tryCatch(utils::read.csv(file.path(d, "tables", "10_MODEL_STATUS.csv"),
                                   stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(st) || !all(c("family", "panel", "status") %in% names(st))) next
    required <- expand.grid(family = r30b_families, panel = r30b_panels,
                            stringsAsFactors = FALSE)
    ok <- vapply(seq_len(nrow(required)), function(i) {
      z <- st[st$family == required$family[i] & st$panel == required$panel[i], , drop = FALSE]
      nrow(z) == 1L && grepl("^PASS", z$status)
    }, logical(1))
    if (all(ok)) return(normalizePath(d, mustWork = TRUE))
  }
  stop("No hay corrida completa con los 9 modelos. Define MT2_30_ROC_SOURCE_DIR para inspeccionar.",
       call. = FALSE)
}

r30b_src <- r30b_choose_source()
r30b_tab <- file.path(r30b_src, "tables")
r30b_out <- file.path(r30b_src, "figures", "ROC_COMPARATIVAS_OVR")
dir.create(r30b_out, recursive = TRUE, showWarnings = FALSE)
r30b_rocs <- utils::read.csv(file.path(r30b_tab, "23_ROC_CURVE_POINTS_OVR.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
r30b_auc <- utils::read.csv(file.path(r30b_tab, "25_ROC_PR_OVR_BY_CLASS.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
r30b_status <- utils::read.csv(file.path(r30b_tab, "10_MODEL_STATUS.csv"),
                                stringsAsFactors = FALSE)
r30b_preds <- utils::read.csv(file.path(r30b_tab, "15_PREDICTIONS_PER_SAMPLE.csv"),
                                 stringsAsFactors = FALSE, check.names = FALSE)
r30b_assert(all(c("family", "panel", "split", "sample_id", "true") %in% names(r30b_preds)),
            "15_PREDICTIONS_PER_SAMPLE.csv no contiene identificadores y clase real")
r30b_preds <- r30b_preds[grepl("^TEST", r30b_preds$split), , drop = FALSE]
for (fam in r30b_families) {
  reference <- NULL
  for (panel in r30b_panels) {
    pred <- r30b_preds[r30b_preds$family == fam & r30b_preds$panel == panel,
                        c("sample_id", "true"), drop = FALSE]
    r30b_assert(nrow(pred) > 0L && !anyNA(pred$sample_id) &&
                  !anyDuplicated(pred$sample_id),
                paste("Predicciones TEST incompletas/duplicadas:", fam, panel))
    pred <- pred[order(pred$sample_id), , drop = FALSE]
    rownames(pred) <- NULL
    if (is.null(reference)) reference <- pred
    else r30b_assert(identical(reference, pred),
             paste("La poblacion TEST difiere entre paneles:", fam, panel))
  }
}
r30b_assert(all(c("family", "panel", "split", "class", "fpr", "tpr") %in% names(r30b_rocs)),
            "23_ROC_CURVE_POINTS_OVR.csv sin columnas ROC esperadas")
r30b_assert(all(c("family", "panel", "split", "class", "roc_auc", "n_pos", "n_neg") %in% names(r30b_auc)),
            "25_ROC_PR_OVR_BY_CLASS.csv sin columnas esperadas")

r30b_rocs <- r30b_rocs[grepl("^TEST", r30b_rocs$split), , drop = FALSE]
r30b_auc <- r30b_auc[grepl("^TEST", r30b_auc$split), , drop = FALSE]
r30b_summary <- list()
r30b_test_splits <- character()
for (fam in r30b_families) for (cls in r30b_classes) {
  sub <- r30b_auc[r30b_auc$family == fam & r30b_auc$class == cls &
                    r30b_auc$panel %in% r30b_panels, , drop = FALSE]
  r30b_assert(nrow(sub) == 3L && setequal(sub$panel, r30b_panels),
               paste("Se requieren 3 AUC para", fam, cls, "(TEST)"))
  r30b_assert(length(unique(sub$n_pos)) == 1L && length(unique(sub$n_neg)) == 1L,
               paste("Distinto número de casos entre modelos:", fam, cls))
  r30b_assert(length(unique(sub$split)) == 1L,
               paste("Splits TEST inconsistentes:", fam, cls))
  r30b_test_splits <- c(r30b_test_splits, unique(sub$split))
  for (panel in r30b_panels) {
    st <- r30b_status[r30b_status$family == fam & r30b_status$panel == panel, , drop = FALSE]
    r30b_assert(nrow(st) == 1L && grepl("^PASS", st$status),
                paste("Modelo no PASS:", fam, panel))
    pts <- r30b_rocs[r30b_rocs$family == fam & r30b_rocs$class == cls &
                       r30b_rocs$panel == panel, , drop = FALSE]
    r30b_assert(nrow(pts) >= 2L && length(unique(pts$split)) == 1L &&
                  identical(unique(pts$split), unique(sub$split)) &&
                  all(is.finite(pts$fpr) & is.finite(pts$tpr)) &&
                  all(pts$fpr >= 0 & pts$fpr <= 1 & pts$tpr >= 0 & pts$tpr <= 1),
                 paste("Curva ROC incompleta:", fam, cls, panel))
  }
  sub <- sub[match(r30b_panels, sub$panel), , drop = FALSE]
  sub$model_label <- unname(r30b_labels[sub$panel])
  r30b_summary[[paste(fam, cls, sep = "_")]] <- sub
}
r30b_summary <- do.call(rbind, r30b_summary)
rownames(r30b_summary) <- NULL
utils::write.csv(r30b_summary,
                 file.path(r30b_tab, "34_ROC_COMPARACION_TRES_MODELOS_OVR_TEST.csv"),
                 row.names = FALSE, na = "", fileEncoding = "UTF-8")

r30b_draw <- function(fam, cls, add_legend = TRUE, size = 1) {
  z <- r30b_summary[r30b_summary$family == fam & r30b_summary$class == cls, , drop = FALSE]
  n_pos <- z$n_pos[[1L]]; n_neg <- z$n_neg[[1L]]
  graphics::plot(NA_real_, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
                 xlab = "Tasa de falsos positivos (1 - especificidad)",
                 ylab = "Sensibilidad", main = paste(fam, "|", cls, "vs resto"),
                 cex.main = 1.1 * size, cex.lab = .92 * size)
  graphics::abline(0, 1, col = "grey60", lty = 3)
  legends <- character()
  for (panel in r30b_panels) {
    d <- r30b_rocs[r30b_rocs$family == fam & r30b_rocs$class == cls &
                      r30b_rocs$panel == panel, , drop = FALSE]
    # Se conserva el orden de umbrales tal como lo exportó el script 30.
    graphics::lines(d$fpr, d$tpr, type = "l", lwd = 3, col = unname(r30b_colors[panel]))
    a <- z$roc_auc[z$panel == panel]
    legends <- c(legends, sprintf("%s  |  AUC %.3f", r30b_labels[panel], a))
  }
  if (add_legend) graphics::legend("bottomright", legend = legends,
    col = unname(r30b_colors[r30b_panels]), lwd = 3, bty = "n", cex = .70 * size,
    bg = "white")
  graphics::mtext(sprintf("TEST reutilizado: %d positivos / %d negativos. Interpretacion exploratoria.",
                          n_pos, n_neg), side = 3, line = .15, cex = .72 * size)
}

# Nueve figuras: por familia (P1, P2, P3) y clase (NP, FD, SO).
for (fam in r30b_families) for (cls in r30b_classes) {
  dest <- file.path(r30b_out, sprintf("ROC_%s_%s_vs_rest_47_FINAL28_MOFA_TEST.png", fam, cls))
  grDevices::png(dest, width = 1950, height = 1550, res = 175, bg = "white")
  tryCatch({
    graphics::par(mar = c(5.4, 5.5, 4.3, 1.5))
    r30b_draw(fam, cls)
  }, finally = grDevices::dev.off())
}

# Tres figuras resumen, una por familia: NP, FD, SO lado a lado.
for (fam in r30b_families) {
  dest <- file.path(r30b_out, sprintf("ROC_%s_NP_FD_SO_TRES_MODELOS_TEST.png", fam))
  grDevices::png(dest, width = 3500, height = 1220, res = 170, bg = "white")
  tryCatch({
    graphics::par(mfrow = c(1, 3), mar = c(5.5, 5.2, 4.0, .8))
    for (cls in r30b_classes) r30b_draw(fam, cls, size = .85)
  }, finally = grDevices::dev.off())
}

# Tres figuras resumen, una por clase: P1, P2 y P3 lado a lado.
for (cls in r30b_classes) {
  dest <- file.path(r30b_out, sprintf("ROC_%s_P1_P2_P3_TRES_MODELOS_TEST.png", cls))
  grDevices::png(dest, width = 3500, height = 1220, res = 170, bg = "white")
  tryCatch({
    graphics::par(mfrow = c(1, 3), mar = c(5.5, 5.2, 4.0, .8))
    for (fam in r30b_families) r30b_draw(fam, cls, size = .85)
  }, finally = grDevices::dev.off())
}

cat("ROC_COMPARATIVAS_SOURCE=", r30b_src, "\n", sep = "")
cat("ROC_COMPARATIVAS_FIGURES=", r30b_out, "\n", sep = "")
cat("ROC_COMPARATIVAS_STATUS=PASS; 9 individuales + 6 resúmenes. No hay ajustes nuevos.\n")

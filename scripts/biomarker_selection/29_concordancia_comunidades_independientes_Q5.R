#!/usr/bin/env Rscript
# ==============================================================================
# 29 - CONCORDANCIA EXPLORATORIA: COMUNIDADES INDEPENDIENTES vs PROGRAMAS Q5
# ==============================================================================
# Script AUTONOMO. No hace source() del paso 28 ni reentrena MOFA.
# Trabaja sobre las tablas CSV exportadas por el paso 28.
# Los programas Q5 son la REFERENCIA CONGELADA, no un resultado del clustering.
#
# Analisis A: comunidades de FEATURES por Spearman OBSERVED y LATENT.
# Analisis B: comunidades de PATHWAYS por similitud de sus perfiles Q5.
# Clustering jerarquico average: K explorado con silhouette y robustez.
# K se selecciona por silhouette SIN mirar la concordancia Q5.
# Comparacion posterior: ARI, NMI, contingencia, composicion y mezclas Q5.
#
# IMPORTANTE: no es bootstrap biologico ni validacion externa. Para features,
# robustez = perturbacion/eliminacion del 15% de las correlaciones conocidas.
# Para pathways, robustez = remuestreo de las filas (features) de la matriz Q5.
# Los programas Q5 y las asociaciones Q5 NO son fuentes independientes.
# Los resultados son DESCRIPTIVOS y generan hipotesis, no las confirman.
# ============================================================================

# ========================= CONFIGURACION ====================================
SELECTION_RUN <- path.expand(
  "~/Documents/PhD/metatest/results_reproducible_20261010/comparacion_modelamientos/upstream_comun/downstream/seleccion_mejor_escenario_train/20261011T233118Z_pid679971"
)

# Dejar vacio para detectar automaticamente la ultima exportacion valida.
# O escribir la ruta de una carpeta concreta '*_density_media' del script 28.
EXPORT_RUN_DIR <- ""

FAMILIES <- c("P1", "P2", "P3")
PATTERNS <- c("FD_like_high", "NP_like_high", "SO_like_high")
K_MAX <- 8L                    # Explorar 2,...,min(K_MAX, n_nodos-1).
N_REPLICATES <- 50L           # Pruebas de robustez: 50 por particion y K.
KEEP_KNOWN_CORR <- 0.85       # Retener 85% de correlaciones en cada replica.
SEED <- 20261020L
PNG_RES <- 180L               # Resolucion de las figuras PNG.
stopifnot(K_MAX >= 2L, N_REPLICATES >= 0L,
          KEEP_KNOWN_CORR > 0, KEEP_KNOWN_CORR <= 1)

# ================ LECTURA / COMPROBACION ====================================
options(stringsAsFactors = FALSE)
set.seed(SEED)
if (!dir.exists(SELECTION_RUN)) {
  stop("No existe SELECTION_RUN: ", SELECTION_RUN, "\nActualiza la ruta del encabezado.")
}

folder_ok <- function(d) {
  file.exists(file.path(d, "00_GLOBAL", "tables", "Q5_FROZEN_ANNOTATION_CATALOG.csv")) &&
    any(file.exists(file.path(d, "P1", PATTERNS,
      "HEATMAPS", "00_SOURCE_TABLES", "tables", "feature_feature_observed_latent_rho.csv")))
}

root <- file.path(SELECTION_RUN, "FINAL_STORY_FIGURES_TABLES")
if (!nzchar(EXPORT_RUN_DIR)) {
  cand <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  cand <- cand[vapply(cand, folder_ok, logical(1))]
  if (!length(cand)) {
    stop("No se encontro una exportacion del script 28 con tablas Q5 y correlaciones en: ", root,
         "\nEjecuta el script 28 hasta su exportacion o establece EXPORT_RUN_DIR.")
  }
  EXPORT_RUN_DIR <- cand[order(file.info(cand)$mtime, decreasing = TRUE)][1L]
}
EXPORT_RUN_DIR <- normalizePath(EXPORT_RUN_DIR, mustWork = TRUE)
if (!folder_ok(EXPORT_RUN_DIR)) {
  stop("La carpeta indicada no tiene tablas Q5 / correlaciones: ", EXPORT_RUN_DIR)
}

OUT_ROOT <- file.path(SELECTION_RUN, "AUDITORIA_CONCORDANCIA_Q5")
OUT_DIR <- file.path(OUT_ROOT, format(Sys.time(), "%Y%m%dT%H%M%S"))
if (dir.exists(OUT_DIR)) OUT_DIR <- paste0(OUT_DIR, "_", Sys.getpid())
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

csv_read <- function(path, required = character()) {
  if (!file.exists(path)) return(NULL)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Columnas ausentes en ", path, ": ", paste(missing, collapse = ", "))
  x
}
csv_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

# ================ MEDIDAS SIN PAQUETES ADICIONALES ===========================
comb2 <- function(z) z * (z - 1) / 2

# Adjusted Rand Index. Invariante a la permutacion de etiquetas.
ari <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- as.character(a[ok]); b <- as.character(b[ok])
  if (length(a) < 3L || length(unique(a)) < 2L || length(unique(b)) < 2L)
    return(NA_real_)
  t <- table(a, b)
  n <- sum(t)
  ab <- sum(comb2(t))
  a2 <- sum(comb2(rowSums(t)))
  b2 <- sum(comb2(colSums(t)))
  expected <- a2 * b2 / comb2(n)
  upper <- (a2 + b2) / 2
  den <- upper - expected
  if (!is.finite(den) || abs(den) < 1e-12) return(NA_real_)
  (ab - expected) / den
}

# Normalized Mutual Information (geometrica), no depende del nombre de cluster.
nmi <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- as.character(a[ok]); b <- as.character(b[ok])
  if (length(a) < 3L || length(unique(a)) < 2L || length(unique(b)) < 2L)
    return(NA_real_)
  p <- prop.table(table(a, b)); pa <- rowSums(p); pb <- colSums(p)
  ind <- which(p > 0, arr.ind = TRUE)
  mi <- sum(p[ind] * log(p[ind] / (pa[ind[, 1L]] * pb[ind[, 2L]])))
  ha <- -sum(pa[pa > 0] * log(pa[pa > 0]))
  hb <- -sum(pb[pb > 0] * log(pb[pb > 0]))
  if (ha <= 0 || hb <= 0) return(NA_real_)
  max(0, min(1, mi / sqrt(ha * hb)))
}

# Silhouette sobre DISTANCIA (no sobre correlacion). Singleton -> silhouette 0.
mean_silhouette <- function(dmat, clusters) {
  n <- length(clusters)
  if (n < 3L || length(unique(clusters)) < 2L || nrow(dmat) != n) return(NA_real_)
  groups <- split(seq_len(n), clusters)
  s <- numeric(n)
  for (idx in groups) {
    if (length(idx) == 1L) next
    others <- groups[!vapply(groups, identical, logical(1), y = idx)]
    for (i in idx) {
      a <- mean(dmat[i, setdiff(idx, i)])
      b <- min(vapply(others, function(j) mean(dmat[i, j]), numeric(1)))
      if (max(a, b) > 0) s[i] <- (b - a) / max(a, b)
    }
  }
  mean(s)
}

# Matriz Rho: solo pares conocidos; los pares no registrados se ponen a rho=0
# EXCLUSIVAMENTE como convencion neutra para construir la distancia.
feature_distances <- function(ids, pairs, col_rho) {
  n <- length(ids)
  d <- matrix(1, n, n, dimnames = list(ids, ids))
  diag(d) <- 0
  known <- matrix(FALSE, n, n)
  i <- match(as.character(pairs$feature1), ids)
  j <- match(as.character(pairs$feature2), ids)
  vals <- suppressWarnings(as.numeric(pairs[[col_rho]]))
  ok <- which(!is.na(i) & !is.na(j) & i != j & is.finite(vals))
  for (h in ok) {
    dd <- 1 - max(-1, min(1, vals[h]))
    d[i[h], j[h]] <- dd; d[j[h], i[h]] <- dd
    known[i[h], j[h]] <- TRUE; known[j[h], i[h]] <- TRUE
  }
  upper <- upper.tri(known)
  list(d = d, known = known,
       coverage = if (sum(upper) > 0) sum(known[upper]) / sum(upper) else 0)
}

# Perfiles Q5: fila=feature, columna=pathway; 0 es ausencia de enlace Q5.
make_association_matrix <- function(links) {
  fids <- sort(unique(as.character(links$feature_id)))
  pids <- paste(as.character(links$library), as.character(links$pathway), sep = "|")
  paths <- sort(unique(pids))
  m <- matrix(0, length(fids), length(paths), dimnames = list(fids, paths))
  for (h in seq_len(nrow(links))) {
    v <- suppressWarnings(as.numeric(links$association_percentile[h]))
    if (is.finite(v)) m[as.character(links$feature_id[h]), pids[h]] <-
      max(m[as.character(links$feature_id[h]), pids[h]], v)
  }
  # Suprimir pathways constantes: la Spearman no estaria definida.
  keep <- apply(m, 2L, function(z) length(unique(z)) > 1L)
  m[, keep, drop = FALSE]
}

# Similitud por perfil entre pathways a partir de la matriz de asociaciones.
pathway_distances <- function(a) {
  s <- suppressWarnings(stats::cor(a, method = "spearman", use = "pairwise.complete.obs"))
  s[!is.finite(s)] <- 0
  diag(s) <- 1
  d <- 1 - s
  d[d < 0] <- 0
  d[d > 2] <- 2
  diag(d) <- 0
  d
}

# Cluster-robustness desde perturbaciones SIN datos biologicos crudos.
# Para features: se omite un subconjunto de correlaciones conocidas.
# Para pathways: se remuestrean FEATURES del perfil de asociaciones Q5.
run_clustering <- function(d, q5, ids, labels, mode,
                           known = NULL, association = NULL) {
  n <- length(ids)
  kvals <- if (n >= 3L) seq.int(2L, min(as.integer(K_MAX), n - 1L)) else integer()
  if (!length(kvals)) return(NULL)
  hc <- stats::hclust(stats::as.dist(d), method = "average")
  base <- lapply(kvals, function(k) stats::cutree(hc, k = k))
  names(base) <- as.character(kvals)
  sil <- vapply(base, function(z) mean_silhouette(d, z), numeric(1))

  stability <- matrix(NA_real_, nrow = N_REPLICATES, ncol = length(kvals))
  if (N_REPLICATES > 0L) for (b in seq_len(N_REPLICATES)) {
    if (mode == "features") {
      db <- d
      ij <- which(known & upper.tri(known), arr.ind = TRUE)
      if (nrow(ij)) {
        nr <- min(nrow(ij), round(nrow(ij) * (1 - KEEP_KNOWN_CORR)))
        drop <- if (nr > 0) sample(seq_len(nrow(ij)), nr) else integer()
        if (length(drop)) {
          ii <- ij[drop, 1L]; jj <- ij[drop, 2L]
          db[cbind(ii, jj)] <- 1
          db[cbind(jj, ii)] <- 1
        }
      }
    } else {
      ind <- sample(seq_len(nrow(association)), size = nrow(association), replace = TRUE)
      db <- pathway_distances(association[ind, , drop = FALSE])
    }
    hc_b <- stats::hclust(stats::as.dist(db), method = "average")
    for (j in seq_along(kvals)) stability[b, j] <- ari(base[[j]], stats::cutree(hc_b, k = kvals[j]))
  }
  stab <- if (N_REPLICATES > 0L) apply(stability, 2L, function(z) {
    if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  }) else rep(NA_real_, length(kvals))

  q5 <- as.character(q5)
  # NA = afiliacion no exclusiva o no disponible: no se fuerza un Q5 artificial.
  ari_q5 <- vapply(base, function(z) ari(z, q5), numeric(1))
  nmi_q5 <- vapply(base, function(z) nmi(z, q5), numeric(1))
  metrics <- data.frame(k = kvals, silhouette = as.numeric(sil),
                        robustness_mean_ARI = as.numeric(stab),
                        ARI_vs_Q5 = as.numeric(ari_q5),
                        NMI_vs_Q5 = as.numeric(nmi_q5),
                        n_nodes = n, n_q5_comparable = sum(!is.na(q5)),
                        n_q5_groups_comparable = length(unique(q5[!is.na(q5)])))
  # Criterio predefinido: silhouette. En caso de empate, K menor.
  candidate <- which(is.finite(metrics$silhouette))
  best <- candidate[which.max(metrics$silhouette[candidate])]
  kbest <- metrics$k[best]
  selected <- base[[as.character(kbest)]]
  assignments <- data.frame(node_id = ids, node_label = as.character(labels),
                            original_Q5_program = q5,
                            independent_community = paste0("C", selected),
                            stringsAsFactors = FALSE)
  # La composicion usa solo comparables; deja patente la mezcla de Q5.
  comp <- data.frame()
  known_ids <- !is.na(q5) & nzchar(q5)
  for (cl in unique(selected)) {
    idx <- which(selected == cl)
    subset_q <- q5[idx][!is.na(q5[idx]) & nzchar(q5[idx])]
    tab <- sort(table(subset_q), decreasing = TRUE)
    if (!length(tab)) {
      comp <- rbind(comp, data.frame(community = paste0("C", cl), n_total = length(idx),
                n_Q5_comparable = 0L, q5_program = NA_character_, n_in_Q5 = 0L,
                fraction_within_community = NA_real_, n_distinct_Q5 = 0L))
    } else {
      for (j in seq_along(tab)) comp <- rbind(comp, data.frame(
        community = paste0("C", cl), n_total = length(idx),
        n_Q5_comparable = length(subset_q), q5_program = names(tab)[j],
        n_in_Q5 = as.integer(tab[j]), fraction_within_community = as.numeric(tab[j] / length(subset_q)),
        n_distinct_Q5 = length(tab)))
    }
  }
  contingency <- if (sum(known_ids)) as.data.frame.matrix(
    table(community = paste0("C", selected[known_ids]), Q5 = q5[known_ids])) else data.frame()
  list(metrics = metrics, selected = assignments, composition = comp,
       contingency = contingency, kbest = kbest, baseline = base)
}

plot_metrics <- function(metrics, dest, label) {
  grDevices::png(dest, width = 2900, height = 1250, res = PNG_RES, bg = "white")
  op <- graphics::par(mfrow = c(1, 3), mar = c(5, 5, 4, 1), oma = c(0, 0, 2, 0))
  on.exit({graphics::par(op); grDevices::dev.off()}, add = TRUE)
  xx <- metrics$k
  graphics::plot(xx, metrics$silhouette, type = "b", lwd = 2, pch = 19,
                 ylim = c(-1, 1), xlab = "K (comunidades)", ylab = "Silhouette media",
                 main = "Criterio no supervisado")
  graphics::abline(h = 0, lty = 3)
  graphics::plot(xx, metrics$robustness_mean_ARI, type = "b", lwd = 2, pch = 19,
                 ylim = c(-1, 1), xlab = "K (comunidades)", ylab = "ARI medio replicas",
                 main = "Robustez (sensibilidad)")
  graphics::abline(h = 0, lty = 3)
  graphics::plot(xx, metrics$ARI_vs_Q5, type = "b", lwd = 2, pch = 19,
                 ylim = c(-1, 1), xlab = "K (comunidades)", ylab = "ARI con Q5",
                 main = "Concordancia posterior Q5")
  graphics::abline(h = 0, lty = 3)
  graphics::mtext(label, side = 3, outer = TRUE, line = 0.4, cex = 0.9)
}

plot_contingency <- function(assignments, dest, label) {
  keep <- !is.na(assignments$original_Q5_program) & nzchar(assignments$original_Q5_program)
  if (!sum(keep)) return(invisible(NULL))
  tab <- table(assignments$original_Q5_program[keep], assignments$independent_community[keep])
  q5_names <- rownames(tab)
  q5_short <- sub("^.*__(P[0-9]+)$", "\\1", q5_names)
  if (anyDuplicated(q5_short)) q5_short <- q5_names
  rownames(tab) <- q5_short
  grDevices::png(dest, width = 2600, height = 1300, res = PNG_RES, bg = "white")
  op <- graphics::par(mar = c(5, 5, 4, 11), xpd = TRUE)
  on.exit({graphics::par(op); grDevices::dev.off()}, add = TRUE)
  prop <- prop.table(tab, margin = 2L)
  colors <- grDevices::hcl.colors(nrow(tab), palette = "Dark 3")
  graphics::barplot(prop, col = colors, border = "white", ylim = c(0, 1),
                    ylab = "Fraccion dentro de cada comunidad", xlab = "Comunidad descubierta",
                    main = label)
  graphics::legend("topright", inset = c(-0.28, 0), legend = rownames(tab),
                   fill = colors, bty = "n", title = "Programa Q5")
}

# ===================== ITERAR FAMILIAS / PATRONES ============================
RESULTS <- list(); STATUS <- list(); ASSIGNMENTS <- list()
add_status <- function(fam, pat, type, signal, state, message_text) {
  STATUS[[length(STATUS) + 1L]] <<- data.frame(
    family = fam, pattern = pat, type = type, signal = signal,
    state = state, detail = message_text, stringsAsFactors = FALSE)
}

for (fam in FAMILIES) for (pat in PATTERNS) {
  prefix <- file.path(EXPORT_RUN_DIR, fam, pat, "HEATMAPS", "00_SOURCE_TABLES", "tables")
  links <- csv_read(file.path(prefix, "feature_pathway_Q5_links.csv"),
                    c("feature_id", "feature_label", "program_id", "library", "pathway", "association_percentile"))
  pairs <- csv_read(file.path(prefix, "feature_feature_observed_latent_rho.csv"),
                    c("feature1", "feature2", "rho_observed", "rho_latent"))
  if (is.null(links) || is.null(pairs) || !nrow(links) || !nrow(pairs)) {
    add_status(fam, pat, "input", "-", "SKIP", "Faltan CSV de correlaciones o links Q5")
    next
  }
  message("Analizando ", fam, " | ", pat)
  links$feature_id <- as.character(links$feature_id)
  links$program_id <- as.character(links$program_id)
  f_ids <- sort(unique(links$feature_id))
  f_labels <- vapply(f_ids, function(fid) {
    z <- unique(links$feature_label[links$feature_id == fid]); as.character(z[1L])
  }, character(1))
  fq <- lapply(f_ids, function(fid) unique(links$program_id[links$feature_id == fid]))
  f_q5 <- vapply(fq, function(z) if (length(z) == 1L && !is.na(z)) z else NA_character_, character(1))
  if (length(unique(f_q5[!is.na(f_q5)])) < 2L) add_status(fam, pat, "features", "Q5", "WARN",
       "Solo un programa Q5 comparable: ARI/NMI vs Q5 no definidos")
  if (sum(is.na(f_q5))) add_status(fam, pat, "features", "Q5", "WARN",
          paste(sum(is.na(f_q5)), "features con mas de un Q5; excluidas solo de ARI/NMI"))

  for (sg in c("observed", "latent")) {
    col_rho <- paste0("rho_", sg)
    obj <- feature_distances(f_ids, pairs, col_rho)
    if (obj$coverage == 0) {
      add_status(fam, pat, "features", sg, "SKIP", "Ninguna correlacion conocida entre las features del escenario")
      next
    }
    if (obj$coverage < 0.70) add_status(fam, pat, "features", sg, "WARN",
         paste("Cobertura correlaciones:", round(100 * obj$coverage, 1), "% (ausentes -> rho=0, solo para clustering)"))
    out <- run_clustering(obj$d, f_q5, f_ids, f_labels, "features", known = obj$known)
    if (is.null(out)) {
      add_status(fam, pat, "features", sg, "SKIP", "Menos de 3 nodos")
      next
    }
    rel <- paste(fam, pat, "features", sg, sep = "_")
    subdir <- file.path(OUT_DIR, fam, pat, paste0("FEATURES_", toupper(sg)))
    dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
    csv_write(out$metrics, file.path(subdir, "01_METRICAS_POR_K.csv"))
    csv_write(out$selected, file.path(subdir, "02_ASIGNACIONES.csv"))
    csv_write(out$composition, file.path(subdir, "03_COMPOSICION_Q5.csv"))
    csv_write(out$contingency, file.path(subdir, "04_CONTINGENCIA.csv"))
    plot_metrics(out$metrics, file.path(subdir, "05_METRICAS.png"), rel)
    plot_contingency(out$selected, file.path(subdir, "06_CONCORDANCIA_Q5.png"), rel)
    chosen <- out$metrics[out$metrics$k == out$kbest, , drop = FALSE]
    chosen$family <- fam; chosen$pattern <- pat
    chosen$type <- "features"; chosen$signal <- sg
    chosen$feature_corr_coverage <- obj$coverage
    RESULTS[[length(RESULTS) + 1L]] <- chosen
    ASSIGNMENTS[[rel]] <- out$selected
    add_status(fam, pat, "features", sg, "OK", paste0("K=", out$kbest))
  }

  # PATHWAYS: una sola matriz de perfiles Q5, independiente del modo OBS/LATENT.
  a <- make_association_matrix(links)
  pids <- colnames(a)
  if (length(pids) < 3L || nrow(a) < 3L) {
    add_status(fam, pat, "pathways", "Q5_PROFILE", "SKIP", "Pocos pathways variables o features")
    next
  }
  pids_raw <- paste(links$library, links$pathway, sep = "|")
  pq <- lapply(pids, function(p) unique(links$program_id[pids_raw == p]))
  p_q5 <- vapply(pq, function(z) if (length(z) == 1L && !is.na(z)) z else NA_character_, character(1))
  if (length(unique(p_q5[!is.na(p_q5)])) < 2L) add_status(fam, pat, "pathways", "Q5_PROFILE", "WARN",
       "Solo un programa Q5 comparable: ARI/NMI vs Q5 no definidos")
  if (sum(is.na(p_q5))) add_status(fam, pat, "pathways", "Q5_PROFILE", "WARN",
       paste(sum(is.na(p_q5)), "pathways afiliados a varios Q5: excluidos solo de ARI/NMI"))
  obj <- run_clustering(pathway_distances(a), p_q5, pids, pids, "pathways", association = a)
  if (!is.null(obj)) {
    rel <- paste(fam, pat, "pathways", "Q5_PROFILE", sep = "_")
    subdir <- file.path(OUT_DIR, fam, pat, "PATHWAYS_Q5_PROFILE")
    dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
    csv_write(obj$metrics, file.path(subdir, "01_METRICAS_POR_K.csv"))
    csv_write(obj$selected, file.path(subdir, "02_ASIGNACIONES.csv"))
    csv_write(obj$composition, file.path(subdir, "03_COMPOSICION_Q5.csv"))
    csv_write(obj$contingency, file.path(subdir, "04_CONTINGENCIA.csv"))
    plot_metrics(obj$metrics, file.path(subdir, "05_METRICAS.png"), rel)
    plot_contingency(obj$selected, file.path(subdir, "06_CONCORDANCIA_Q5.png"), rel)
    chosen <- obj$metrics[obj$metrics$k == obj$kbest, , drop = FALSE]
    chosen$family <- fam; chosen$pattern <- pat
    chosen$type <- "pathways"; chosen$signal <- "Q5_PROFILE"
    chosen$feature_corr_coverage <- NA_real_
    RESULTS[[length(RESULTS) + 1L]] <- chosen
    add_status(fam, pat, "pathways", "Q5_PROFILE", "OK", paste0("K=", obj$kbest))
  }
}

if (!length(RESULTS)) {
  csv_write(do.call(rbind, STATUS), file.path(OUT_DIR, "00_ESTADO.csv"))
  stop("No hubo combinaciones evaluables. Revisa 00_ESTADO.csv en ", OUT_DIR)
}
final <- do.call(rbind, RESULTS)
final <- final[, c("family", "pattern", "type", "signal", "k", "n_nodes", "n_q5_comparable",
                   "n_q5_groups_comparable", "silhouette", "robustness_mean_ARI",
                   "ARI_vs_Q5", "NMI_vs_Q5", "feature_corr_coverage")]
csv_write(final, file.path(OUT_DIR, "00_RESUMEN_CONCORDANCIA.csv"))
csv_write(do.call(rbind, STATUS), file.path(OUT_DIR, "00_ESTADO.csv"))

# Concordancia entre comunidades de features OBSERVED y LATENT, sobre mismos nodos.
cross <- list()
for (fam in FAMILIES) for (pat in PATTERNS) {
  key0 <- paste(fam, pat, "features", "observed", sep = "_")
  key1 <- paste(fam, pat, "features", "latent", sep = "_")
  if (is.null(ASSIGNMENTS[[key0]]) || is.null(ASSIGNMENTS[[key1]])) next
  x <- ASSIGNMENTS[[key0]]; y <- ASSIGNMENTS[[key1]]
  ids <- intersect(x$node_id, y$node_id)
  v0 <- x$independent_community[match(ids, x$node_id)]
  v1 <- y$independent_community[match(ids, y$node_id)]
  cross[[length(cross) + 1L]] <- data.frame(family = fam, pattern = pat,
       n_common_features = length(ids), ARI_observed_vs_latent = ari(v0, v1),
       NMI_observed_vs_latent = nmi(v0, v1))
}
if (length(cross)) csv_write(do.call(rbind, cross),
  file.path(OUT_DIR, "00_CONCORDANCIA_OBSERVED_LATENT.csv"))

notes <- c(
  "AUDITORIA EXPLORATORIA DE CONCORDANCIA CON PROGRAMAS Q5",
  paste("Generado:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Fuente exacta:", EXPORT_RUN_DIR),
  "",
  "QUE SE CALCULA",
  "- Features: clustering average sobre 1-rho Spearman, por separado OBSERVED/LATENT.",
  "- Pathways: clustering average sobre 1-rho Spearman de perfiles de asociacion Q5.",
  "- K candidato: 2 hasta min(K_MAX, n_nodos - 1). Se elige por silhouette maxima.",
  "- La concordancia ARI/NMI con Q5 se calcula DESPUES de elegir K; NO interviene en K.",
  "- Composicion: fraccion de cada comunidad perteneciente a uno o varios Q5.",
  "- OBSERVED vs LATENT: ARI de sus comunidades seleccionadas para mismas features.",
  "",
  "LIMITACIONES IMPORTANTES",
  "- Las matrices de correlacion son GLOBAL descriptivas, no TEST externo.",
  "- No se dispone de perfiles de muestras para bootstrap biologico.",
  "- Para features, robustez = eliminar aleatoriamente 15% de pares de rho conocidos; pares ausentes -> rho neutro.",
  "- Para pathways, robustez = bootstrap de features de la matriz de asociacion Q5 (no de muestras).",
  "- Los pathway profiles derivan de Q5: una concordancia alta NO es validacion independiente de Q5.",
  "- El original Q5 se obtuvo de comunidades biologicas previas. Sus identificadores son locales a familia/lado.",
  "- Un pathway o feature con Q5 no exclusivo se excluye SOLO del ARI/NMI, no del clustering.",
  "- Si solo existe un programa Q5, ARI/NMI con Q5 no se pueden interpretar (NA).",
  "- K y silhouette pueden ser poco fiables con pocos nodos, pares rho ausentes o perfiles no independientes.",
  "- Esta auditoria NO investiga comunidades entre FD/NP/SO ni homologias P1/P2/P3: cada patron se analiza solo.",
  "- Las comunidades descubiertas NO sustituyen la anotacion/programas Q5 del script 28.",
  "",
  "SALIDAS",
  "00_RESUMEN_CONCORDANCIA.csv: K seleccionado y ARI/NMI por familia/patron/tipo/senal.",
  "00_CONCORDANCIA_OBSERVED_LATENT.csv: concordancia entre particiones de features.",
  "00_ESTADO.csv: advertencias y combinaciones no evaluables.",
  "En cada carpeta: metricas por K, asignaciones, composicion, contingencia y PNG.",
  "",
  "LECTURA DE RESULTADOS",
  "ARI cercano a 1: mismas particiones salvo etiquetas; cercano a 0: concordancia semejante al azar.",
  "ARI negativo: concordancia inferior a la esperada por azar bajo el modelo nulo de ARI.",
  "Silhouette: compactacion/separacion para la distancia usada; NO equivale a interpretacion biologica.",
  "Robustez: ARI con la particion original tras perturbacion; no intervalo de confianza biologico."
)
writeLines(notes, file.path(OUT_DIR, "LEER_PRIMERO.txt"), useBytes = TRUE)
writeLines(capture.output(utils::sessionInfo()), file.path(OUT_DIR, "SESSION_INFO.txt"))
message("\nAUDITORIA TERMINADA: ", OUT_DIR)
message("Combinaciones evaluadas: ", nrow(final), ". Lee 00_RESUMEN_CONCORDANCIA.csv y LEER_PRIMERO.txt")

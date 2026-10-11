#!/usr/bin/env Rscript
# METATEST: seleccion downstream de modelos MOFA finales, exclusivamente TRAIN.
# Version 2.0.1 | Seleccion por meseta de calidad TRAIN + cobertura biologica multivista.
# R >= 4.1. Solo requiere paquetes base/recommended de R.
# Si MOFA2, vegan o cluster estan instalados, utiliza/verifica sus APIs.
# Ver METODO_Y_USO.md para definiciones, limitaciones y verificacion realizada.

# ======================== CONFIGURACION EDITABLE ========================
UPSTREAM_DIR <- Sys.getenv("UPSTREAM_DIR", unset = paste0(
  "~/Documents/PhD/metatest/results_upstream_local_20260908/",
  "comparacion_modelamientos/upstream_comun"))
DOWNSTREAM_DIR <- Sys.getenv("DOWNSTREAM_DIR", unset = "")
# Vacio: <upstream_comun>/downstream/seleccion_mejor_escenario_train/
FAMILIES <- c("P1", "P2", "P3")  # Tres analisis independientes por defecto.
SEED <- 20260910L
N_PERM <- 9999L  # Global: Monte Carlo. Pairwise: exacto si cabe en este limite.
TOP_N <- 15L
TOP_FACTOR_PLOTS <- 5L
MAKE_PLOTS <- TRUE
Z_GEOMETRY <- "unit_variance"  # "unit_variance" o "native"; todos los factores.
TRAIN_MANIFEST_FILE <- Sys.getenv("TRAIN_MANIFEST_FILE", unset = "")
# Si esta vacio, se utiliza la lista cerrada verificada en el RDS adjunto.
# Un manifiesto alternativo debe contener exactamente columnas sample y grupo.
TRAIN_MANIFEST <- data.frame(
  sample = c("NP_01", "NP_04", "NP_06", "NP_10", "NP_11", "NP_12",
             "FD_25", "FD_28", "FD_29", "SO_13", "SO_15", "SO_16",
             "SO_17", "SO_18", "SO_21", "SO_22"),
  grupo = c(rep("NP", 6), rep("FD", 3), rep("SO", 7)),
  stringsAsFactors = FALSE)
EXPECTED_COUNTS <- c(NP = 6L, FD = 3L, SO = 7L)
VIEWS <- c("transcriptomica", "proteomica", "metabolomica", "clinical")
GROUPS <- c("NP", "FD", "SO")
# MOFA2 guarda variance_explained en porcentaje. No inferir unidades por magnitud.
MOFA_R2_DIVISOR <- 100
# Diferencias absolutas sobre metricas en escala de proporcion/silhouette.
# Margenes operativos predeclarados, NO umbrales de equivalencia estadistica.
QUALITY_COMPONENTS <- c(
  "pairwise_separation_worst",
  "silhouette_worst",
  "silhouette_macro",
  "dispersion_balance",
  "view_balance_entropy"
)

QUALITY_EPS <- 1e-12
NUMERICAL_TOL <- 1e-12
FAMILY_DIRS <- c(P1 = "escenarios", P2 = "escenarios_con_FLI",
                 P3 = "escenarios_SIN_FLI")
FAMILY_PATTERNS <- c(P1 = "^P1_(TODAS|ALL)_qtx[0-9]+_qpr[0-9]+_qme[0-9]+_qcl[0-9]+$",
  P2 = "^P2_CON_FLI_qtx[0-9]+_qpr[0-9]+_qme[0-9]+_qcl[0-9]+$",
  P3 = "^P3_SIN_FLI_qtx[0-9]+_qpr[0-9]+_qme[0-9]+_qcl[0-9]+$")
# ====================== FIN DE CONFIGURACION ============================

assert <- function(ok, msg) { if (!isTRUE(ok)) stop(msg, call. = FALSE) }
bind_rows <- function(xs) {
  xs <- Filter(function(x) !is.null(x) && NROW(x) > 0L, xs)
  if (!length(xs)) return(data.frame())
  xs <- lapply(xs, as.data.frame, stringsAsFactors = FALSE)
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  do.call(rbind, lapply(xs, function(x) {
    for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  }))
}
write_csv <- function(x, path) {
  if (is.null(x) || !ncol(x)) x <- data.frame(status = character())
  utils::write.table(x, path, sep = ",", row.names = FALSE, col.names = TRUE,
                     na = "NA", qmethod = "double", fileEncoding = "UTF-8")
}
finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
clamp01 <- function(x) pmin(1, pmax(0, x))
fmt <- function(x) ifelse(is.finite(x), ifelse(abs(x) > 0 & abs(x) < 1e-4,
  formatC(x, digits = 3, format = "e"), formatC(x, digits = 5, format = "f")), "NA")
pkg_version <- function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not_installed"
}
canonical_path <- function(p) {
  p <- path.expand(p)
  if (file.exists(p) || dir.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
  parent <- dirname(p)
  assert(!identical(parent, p), paste("Ruta no resoluble:", p))
  file.path(canonical_path(parent), basename(p))
}
is_inside <- function(p, root) identical(p, root) || startsWith(p, paste0(root, "/"))
entropy_stats <- function(x, denominator = length(x)) {
  x <- pmax(as.numeric(x), 0)
  assert(all(is.finite(x)), "Contribuciones no finitas para entropia")
  if (sum(x) <= 0) return(c(entropy = 0, effective = 0, max_share = 0))
  p <- x / sum(x)
  h <- -sum(p[p > 0] * log(p[p > 0]))
  c(entropy = if (denominator > 1) h / log(denominator) else 0,
    effective = exp(h), max_share = max(p))
}
get_slot <- function(object, name) {
  # Leer atributos S4 no requiere inventar una definicion local de clase MOFA.
  assert(isS4(object), "El objeto no es S4")
  assert(name %in% names(attributes(object)), paste("Falta slot", name))
  attr(object, name, exact = TRUE)
}

validate_manifest <- function(manifest) {
  assert(all(c("sample", "grupo") %in% names(manifest)), "Manifiesto sin sample/grupo")
  manifest <- manifest[, c("sample", "grupo"), drop = FALSE]
  manifest[] <- lapply(manifest, as.character)
  assert(!anyNA(manifest) && all(nzchar(manifest$sample)), "Manifiesto con valores vacios")
  assert(!anyDuplicated(manifest$sample), "Manifiesto con IDs duplicados")
  assert(setequal(unique(manifest$grupo), GROUPS), "Grupos inesperados en manifiesto")
  counts <- table(factor(manifest$grupo, levels = GROUPS))
  assert(identical(as.integer(counts), as.integer(EXPECTED_COUNTS[GROUPS])),
         "TRAIN debe contener NP=6, FD=3, SO=7; adapte EXPECTED_COUNTS solo si cambia el estudio")
  manifest
}

audit_data_samples <- function(object, sample_ids, factor_names) {
  dat <- get_slot(object, "data")
  assert(is.list(dat) && setequal(names(dat), VIEWS), "Las vistas de @data no son las cuatro esperadas")
  w <- get_slot(object, "expectations")$W
  assert(is.list(w) && setequal(names(w), VIEWS), "Las vistas de W no coinciden")
  feature_n <- setNames(integer(length(VIEWS)), VIEWS)
  data_kind <- character()
  for (v in VIEWS) {
    assert(is.list(dat[[v]]) && identical(names(dat[[v]]), "group1"),
           paste(v, "debe contener un unico grupo MOFA group1"))
    x <- dat[[v]][[1L]]
    if ("matrix_placeholder" %in% attr(x, "class", exact = TRUE)) {
      fetch <- function(nm) if (isS4(x)) attr(x, nm, exact = TRUE) else x[[nm]]
      ids <- as.character(fetch("colnames"))
      nr <- as.integer(fetch("nrow")); nc <- as.integer(fetch("ncol"))
      features <- as.character(fetch("rownames"))
      data_kind <- c(data_kind, "matrix_placeholder")
    } else {
      assert(length(dim(x)) == 2L, paste("Matriz no reconocida:", v))
      ids <- colnames(x); features <- rownames(x); nr <- nrow(x); nc <- ncol(x)
      assert(all(is.finite(x)), paste("Datos incompletos/no finitos:", v))
      data_kind <- c(data_kind, "matrix")
    }
    assert(length(nr) == 1L && nr > 0L && length(nc) == 1L && nc == length(sample_ids),
           paste("Dimensiones de datos incompatibles con TRAIN:", v))
    assert(length(ids) == length(sample_ids) && !anyDuplicated(ids) && setequal(ids, sample_ids),
           paste("Muestras adicionales/ausentes en datos de", v, "; modelo rechazado completo"))
    assert(is.matrix(w[[v]]) && nrow(w[[v]]) == nr &&
      ncol(w[[v]]) == length(factor_names) && all(is.finite(w[[v]])), paste("W invalido:", v))
    assert(setequal(colnames(w[[v]]), factor_names), paste("Factores de W no coinciden:", v))
    assert(length(features) == nr && !anyDuplicated(features) &&
      identical(rownames(w[[v]]), features), paste("Features de W/datos no coinciden:", v))
    feature_n[v] <- nr
  }
  list(feature_n = feature_n, data_storage = paste(unique(data_kind), collapse = ";"))
}

extract_model <- function(path, manifest) {
  object <- readRDS(path)
  assert(isS4(object) && "MOFA" %in% class(object), "El RDS no es S4 de clase MOFA")
  md <- as.data.frame(get_slot(object, "samples_metadata"), stringsAsFactors = FALSE)
  assert(all(c("sample", "grupo", "group") %in% names(md)), "Metadata sin sample/grupo/group")
  md <- md[, c("sample", "grupo", "group"), drop = FALSE]
  md[] <- lapply(md, as.character)
  assert(!anyNA(md) && !anyDuplicated(md$sample), "Metadata con NA o muestra duplicada")
  # No se acepta filtrar TRAIN despues de cargar un modelo con otras muestras.
  assert(nrow(md) == nrow(manifest) && setequal(md$sample, manifest$sample),
         "Muestras distintas de la lista cerrada TRAIN; modelo rechazado sin calcular metricas")
  assert(identical(unique(md$group), "group1"), "Solo se admite un grupo tecnico MOFA group1")
  ix <- match(manifest$sample, md$sample)
  assert(identical(md$grupo[ix], manifest$grupo), "Etiquetas grupo distintas del manifiesto TRAIN")
  md <- md[ix, , drop = FALSE]
  rownames(md) <- NULL
  has_mofa <- requireNamespace("MOFA2", quietly = TRUE)
  zz <- if (has_mofa) MOFA2::get_factors(object, scale = FALSE) else get_slot(object, "expectations")$Z
  assert(is.list(zz) && identical(names(zz), "group1"), "Z no tiene exclusivamente group1")
  z <- zz[[1L]]
  assert(is.matrix(z) && nrow(z) == nrow(manifest) && ncol(z) >= 1L, "Dimensiones de Z invalidas")
  assert(!is.null(rownames(z)) && !anyDuplicated(rownames(z)) &&
    setequal(rownames(z), manifest$sample), "IDs de Z distintos de TRAIN")
  assert(!is.null(colnames(z)) && !anyDuplicated(colnames(z)), "Factores sin nombres unicos")
  z_was_reordered <- !identical(rownames(z), manifest$sample)
  z <- z[manifest$sample, , drop = FALSE]
  assert(all(is.finite(z)), "Z tiene NA/Inf")
  da <- audit_data_samples(object, manifest$sample, colnames(z))
  # En el RDS real, las matrices son placeholders: nunca cargar datos externos.
  cache <- get_slot(object, "cache")
  assert(is.list(cache) && "variance_explained" %in% names(cache),
         "Falta variance_explained guardada; no se reconstruyen datos ni se accede a rutas externas")
  ve <- if (has_mofa) MOFA2::get_variance_explained(object, as.data.frame = FALSE) else cache$variance_explained
  assert(all(c("r2_total", "r2_per_factor") %in% names(ve)), "Cache R2 incompleta")
  assert(identical(names(ve$r2_total), "group1") && identical(names(ve$r2_per_factor), "group1"),
         "Cache R2 no corresponde exclusivamente a group1")
  total <- ve$r2_total[[1L]]
  assert(is.numeric(total) && setequal(names(total), VIEWS), "R2 total sin nombres de vistas validos")
  total <- total[VIEWS] / MOFA_R2_DIVISOR
  pf <- ve$r2_per_factor[[1L]]
  assert(is.matrix(pf) && setequal(rownames(pf), colnames(z)) &&
           setequal(colnames(pf), VIEWS), "R2 factor x vista mal etiquetado")
  pf <- pf[colnames(z), VIEWS, drop = FALSE] / MOFA_R2_DIVISOR
  assert(all(is.finite(total)) && all(is.finite(pf)), "Variance explained no finita")
  assert(all(total <= 1 + 1e-6) && all(pf <= 1 + 1e-6), "R2 mayor que 100%; compruebe unidades/cache")
  centers <- colMeans(z); scales <- apply(z, 2, stats::sd)
  zero_factor <- scales <= NUMERICAL_TOL
  scales[zero_factor] <- 1
  z_eval <- sweep(z, 2, centers, "-")
  if (Z_GEOMETRY == "unit_variance") z_eval <- sweep(z_eval, 2, scales, "/")
  if (any(zero_factor)) z_eval[, zero_factor] <- 0
  assert(sum(z_eval^2) > NUMERICAL_TOL, "Z completamente degenerada")
  list(Z = z_eval, Z_native = z, metadata = md, view_total = total, factor_view = pf,
       center = centers, scale = scales, zero_factor = zero_factor,
       feature_n = da$feature_n, data_storage = da$data_storage,
       backend = if (has_mofa) "MOFA2_API_cached" else "S4_slots_cached_no_MOFA2",
       reordered = z_was_reordered, source = normalizePath(path, mustWork = TRUE),
       source_md5 = unname(tools::md5sum(path)))
}


invert_permutation_matrix <- function(pmat) {
  n <- ncol(pmat)
  inv <- matrix(NA_integer_, nrow = nrow(pmat), ncol = n)
  
  for (i in seq_len(nrow(pmat))) {
    inv[i, pmat[i, ]] <- seq_len(n)
  }
  
  inv
}
# Permutaciones libres conservan los tamanos NP/FD/SO; NO bloquear por grupo.
# La evaluacion balanceada no significa permutar solamente dentro de cada grupo.
make_banks <- function(manifest, n_perm, seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  
  g <- manifest$grupo
  n <- length(g)
  
  # Banco definido conceptualmente como permutacion de ETIQUETAS
  label_idx <- t(replicate(n_perm, sample.int(n)))
  
  # La permutacion equivalente de la RESPUESTA es la inversa
  response_idx <- invert_permutation_matrix(label_idx)
  
  memberships <- lapply(
    GROUPS,
    function(gr)
      matrix(as.numeric(g[label_idx] == gr), nrow = n_perm)
  )
  names(memberships) <- GROUPS
  
  pairs <- list(c("NP", "FD"), c("NP", "SO"), c("FD", "SO"))
  
  pb <- lapply(pairs, function(p) {
    rows <- which(g %in% p)
    gg <- g[rows]
    nn <- length(rows)
    na <- sum(gg == p[1])
    
    n_labels <- choose(nn, na)
    exact <- n_labels <= n_perm + 1L
    
    combinations <- if (exact) {
      utils::combn(nn, na)
    } else {
      replicate(n_perm, sample.int(nn, na))
    }
    
    A <- matrix(0, nrow = ncol(combinations), ncol = nn)
    for (j in seq_len(nrow(A))) {
      A[j, combinations[, j]] <- 1
    }
    
    list(
      groups = p,
      rows = rows,
      A = A,
      exact = exact,
      n_labels = n_labels,
      n_draws = nrow(A),
      n_a = na,
      n_b = nn - na
    )
  })
  
  names(pb) <- vapply(
    pairs, paste, collapse = "_",
    FUN.VALUE = character(1)
  )
  
  list(
    global_label_indices = label_idx,
    global_response_indices = response_idx,
    global_memberships = memberships,
    pairwise = pb,
    sample_order = manifest$sample,
    seed = seed,
    n_perm = n_perm
  )
}
permutation_p <- function(null, observed, exact = FALSE) {
  if (!is.finite(observed) || any(!is.finite(null))) return(NA_real_)
  extreme <- sum(null >= observed - NUMERICAL_TOL)
  if (exact) extreme / length(null) else (extreme + 1) / (length(null) + 1)
}

global_factor_metrics <- function(z, g, banks) {
  zc <- sweep(z, 2, colMeans(z), "-")
  sst <- colSums(zc^2)
  ssb <- numeric(ncol(z))
  null_ssb <- matrix(
    0,
    nrow = nrow(banks$global_label_indices),
    ncol = ncol(z)
  )
  for (gr in GROUPS) {
    ng <- sum(g == gr)
    ssb <- ssb + colSums(zc[g == gr, , drop = FALSE])^2 / ng
    null_ssb <- null_ssb + (banks$global_memberships[[gr]] %*% zc)^2 / ng
  }
  eta <- ifelse(sst > NUMERICAL_TOL, clamp01(ssb / pmax(sst, NUMERICAL_TOL)), 0)
  null_eta <- sweep(null_ssb, 2, pmax(sst, NUMERICAL_TOL), "/")
  p <- vapply(seq_along(eta), function(k) {
    if (sst[k] <= NUMERICAL_TOL) 1 else permutation_p(null_eta[, k], eta[k])
  }, numeric(1))
  fac <- data.frame(factor = colnames(z), r2_group = eta, eta2_group = eta,
    between_group_SS = ssb, total_SS = sst, p_grupo = p,
    q_grupo = stats::p.adjust(p, "BH"), constant_factor = sst <= NUMERICAL_TOL,
    p_method = "free_label_permutation_MC", stringsAsFactors = FALSE)
  global_r2 <- clamp01(sum(ssb) / sum(sst))
  global_p <- permutation_p(rowSums(null_ssb) / sum(sst), global_r2)
  # Verificacion independiente con vegan cuando esta disponible.
  if (requireNamespace("vegan", quietly = TRUE)) {
    d <- stats::dist(z)
    a <- vegan::adonis2(d ~ grupo, data = data.frame(grupo = factor(g)),
                        permutations = banks$global_response_indices,
                        parallel = 1)
    assert(abs(a$R2[1] - global_r2) < 1e-8, "R2 global no coincide con vegan::adonis2")
    assert(
      abs(a$`Pr(>F)`[1] - global_p) <=
        1 / (nrow(banks$global_label_indices) + 1) + 1e-8,
      "p global no coincide con vegan::adonis2"
    )
  }
  list(factors = fac, R2 = global_r2, p = global_p, ssb = ssb)
}

# R2 equilibrado de centroides: B/(B+W), B=||mu_A-mu_B||^2/4,
# W=(tr(S_A)+tr(S_B))/2; S_g usa denominador n_g-1.
# Es un efecto descriptivo de prevalencias iguales, NO el R2 de adonis2.
pair_statistics <- function(z, A) {
  n <- nrow(z); na <- rowSums(A); nb <- n - na
  sums_a <- A %*% z
  sums_b <- matrix(colSums(z), nrow(A), ncol(z), byrow = TRUE) - sums_a
  means_a <- sums_a / na; means_b <- sums_b / nb
  delta2 <- rowSums((means_a - means_b)^2)
  sq <- rowSums(z^2)
  ssa <- pmax(0, as.vector(A %*% sq) - rowSums(sums_a^2) / na)
  ssb <- pmax(0, sum(sq) - as.vector(A %*% sq) - rowSums(sums_b^2) / nb)
  within <- (ssa / (na - 1) + ssb / (nb - 1)) / 2
  between <- delta2 / 4
  separation <- ifelse(between + within > NUMERICAL_TOL, between / (between + within), 0)
  ss_between <- na * nb / n * delta2
  sst <- ss_between + ssa + ssb
  ordinary <- ifelse(sst > NUMERICAL_TOL, ss_between / sst, 0)
  data.frame(separation = clamp01(separation), permanova_R2 = clamp01(ordinary),
             centroid_distance = sqrt(delta2), within_trace_a = ssa / (na - 1),
             within_trace_b = ssb / (nb - 1), between_equal_prior = between,
             within_equal_prior = within)
}
pairwise_metrics <- function(z, g, banks) {
  bind_rows(lapply(names(banks$pairwise), function(nm) {
    b <- banks$pairwise[[nm]]; zz <- z[b$rows, , drop = FALSE]; gg <- g[b$rows]
    obs <- pair_statistics(zz, matrix(as.numeric(gg == b$groups[1]), nrow = 1))
    null <- pair_statistics(zz, b$A)
    data.frame(contrast = nm, group_a = b$groups[1], group_b = b$groups[2],
      n_a = b$n_a, n_b = b$n_b, obs,
      p = permutation_p(null$separation, obs$separation, b$exact),
      permanova_p = permutation_p(null$permanova_R2, obs$permanova_R2, b$exact),
      permutation_method = if (b$exact) "exact_all_unique_label_allocations" else "MC_label_allocations",
      n_permutations = b$n_draws, n_possible_label_allocations = b$n_labels,
      p_resolution = if (b$exact) 1 / b$n_draws else 1 / (b$n_draws + 1),
      stringsAsFactors = FALSE)
  }))
}

silhouette_metrics <- function(z, g, sample) {
  d <- as.matrix(stats::dist(z))
  widths <- vapply(seq_len(nrow(z)), function(i) {
    own <- which(g == g[i] & seq_along(g) != i)
    a <- mean(d[i, own])
    other_means <- vapply(setdiff(GROUPS, g[i]), function(gr) mean(d[i, g == gr]), numeric(1))
    b <- min(other_means)
    if (max(a, b) > 0) (b - a) / max(a, b) else 0
  }, numeric(1))
  if (requireNamespace("cluster", quietly = TRUE)) {
    check <- cluster::silhouette(as.integer(factor(g, levels = GROUPS)), stats::as.dist(d))
    assert(max(abs(check[, "sil_width"] - widths)) < 1e-10, "Silhouette no coincide con cluster")
  }
  by_group <- setNames(vapply(GROUPS, function(gr) mean(widths[g == gr]), numeric(1)), GROUPS)
  list(global = mean(widths), by_group = by_group, macro = mean(by_group), worst = min(by_group),
       samples = data.frame(sample = sample, grupo = g, silhouette = widths))
}

dispersion_metrics <- function(z, g, banks) {
  dd <- bind_rows(lapply(GROUPS, function(gr) {
    x <- z[g == gr, , drop = FALSE]
    radius <- sqrt(rowSums(sweep(x, 2, colMeans(x), "-")^2))
    data.frame(grupo = gr, n = nrow(x), mean_centroid_distance = mean(radius),
      mean_centroid_distance_bias_adjusted = mean(radius) * sqrt(nrow(x) / (nrow(x) - 1)),
      within_trace_unbiased = sum(radius^2) / (nrow(x) - 1),
      within_RMS_bias_adjusted = sqrt(sum(radius^2) / (nrow(x) - 1)))
  }))
  p <- NA_real_; method <- "centroid_radii_only"
  note <- "vegan no disponible: se guardan dispersiones; p_dispersion_global queda NA diagnostico"
  if (requireNamespace("vegan", quietly = TRUE)) {
    note <- ""
    attempt <- tryCatch({
      bd <- vegan::betadisper(stats::dist(z), factor(g, levels = GROUPS),
                             type = "median", bias.adjust = TRUE)
      pp <- vegan::permutest(
        bd,
        permutations = banks$global_response_indices
      )
      as.numeric(pp$tab[1, "Pr(>F)"])
    }, error = function(e) { note <<- conditionMessage(e); NA_real_ })
    p <- attempt; method <- "vegan_betadisper_spatial_median_bias_adjusted"
  }
  list(table = dd, p = p, method = method, note = note)
}
quality_score <- function(x) {
  n <- nrow(x)
  assert(n >= 1L, "quality_score requiere al menos un escenario")

  q <- matrix(
    NA_real_,
    nrow = n,
    ncol = length(QUALITY_COMPONENTS),
    dimnames = list(NULL, QUALITY_COMPONENTS)
  )

  q[, "pairwise_separation_worst"] <- clamp01(as.numeric(x$pairwise_separation_worst))
  q[, "silhouette_worst"] <- clamp01(as.numeric(x$silhouette_worst))
  q[, "silhouette_macro"] <- clamp01(as.numeric(x$silhouette_macro))
  q[, "dispersion_balance"] <- clamp01(as.numeric(x$dispersion_balance))
  q[, "view_balance_entropy"] <- clamp01(as.numeric(x$view_balance_entropy))

  assert(
    is.matrix(q) &&
      identical(dim(q), c(n, length(QUALITY_COMPONENTS))) &&
      all(is.finite(q)),
    "Matriz de componentes de quality_score invalida"
  )

  exp(rowMeans(log(pmax(q, QUALITY_EPS))))
}

detect_quality_plateau <- function(score) {
  
  assert(
    length(score) > 0 &&
      all(is.finite(score)),
    "quality_score invalido"
  )
  
  ord <- order(-score, method = "radix")
  s <- score[ord]
  
  in_plateau <- rep(FALSE, length(score))
  
  if (
    length(s) == 1L ||
    max(s) - min(s) <= NUMERICAL_TOL
  ) {
    in_plateau[] <- TRUE
    
    return(list(
      in_plateau = in_plateau,
      cutoff = min(s),
      largest_drop = 0,
      split_after = length(s)
    ))
  }
  
  drops <- s[-length(s)] - s[-1L]
  
  k <- which.max(drops)
  
  in_plateau[ord[seq_len(k)]] <- TRUE
  
  list(
    in_plateau = in_plateau,
    cutoff = s[k],
    next_score = s[k + 1L],
    largest_drop = drops[k],
    split_after = k
  )
}

biological_coverage <- function(x) {
  feature_cols <- paste0("n_features_", VIEWS)
  n <- nrow(x)
  assert(n >= 1L, "biological_coverage requiere al menos un escenario")

  baseline_idx <- which(x$baseline)
  if (length(baseline_idx) == 1L) {
    reference <- setNames(
      as.numeric(unlist(
        x[baseline_idx, feature_cols, drop = FALSE],
        use.names = FALSE
      )),
      VIEWS
    )
    reference_source <- "q00_baseline"
  } else {
    reference <- setNames(
      vapply(
        feature_cols,
        function(nm) max(as.numeric(x[[nm]]), na.rm = TRUE),
        numeric(1)
      ),
      VIEWS
    )
    reference_source <- "family_maximum_fallback"
  }

  assert(
    length(reference) == length(VIEWS) &&
      all(is.finite(reference)) &&
      all(reference > 0),
    "Referencia de cobertura biologica invalida"
  )

  ratios <- matrix(
    NA_real_,
    nrow = n,
    ncol = length(VIEWS),
    dimnames = list(NULL, VIEWS)
  )

  for (j in seq_along(feature_cols)) {
    ratios[, j] <- as.numeric(x[[feature_cols[j]]]) / reference[j]
  }

  assert(
    is.matrix(ratios) &&
      identical(dim(ratios), c(n, length(VIEWS))) &&
      all(is.finite(ratios)),
    "Matriz de retenciones biologicas invalida"
  )

  assert(
    all(ratios <= 1 + 1e-8),
    "Hay un escenario con mas features que la referencia biologica; revisar q00/upstream"
  )

  # Mantener la matriz como primer argumento de pmin/pmax conserva dim/dimnames.
  ratios <- pmin(ratios, 1)
  ratios <- pmax(ratios, QUALITY_EPS)

  gmean <- exp(rowMeans(log(ratios)))
  minimum <- apply(ratios, 1L, min)

  list(
    score = gmean,
    minimum = minimum,
    ratios = ratios,
    reference = reference,
    reference_source = reference_source
  )
}

view_metrics <- function(m) {
  total <- m$view_total; pf <- m$factor_view
  # Solo las participaciones usadas en entropias se truncaran en cero.
  # Los R2 negativos originales se conservan, se cuentan y se auditan.
  positive <- pmax(total, 0)
  vh <- entropy_stats(positive, 4)
  fh <- t(apply(pf, 1, entropy_stats, denominator = 4))
  factors <- data.frame(factor = rownames(pf), entropy_factor = fh[, "entropy"],
    effective_views_factor = fh[, "effective"], max_view_share = fh[, "max_share"],
    sum_positive_factor_R2 = rowSums(pmax(pf, 0)),
    zero_positive_contribution = rowSums(pmax(pf, 0)) <= 0, stringsAsFactors = FALSE)
  fv <- bind_rows(lapply(seq_len(nrow(pf)), function(k) {
    p <- pmax(pf[k, ], 0); if (sum(p) > 0) p <- p / sum(p)
    data.frame(factor = rownames(pf)[k], view = VIEWS, R2 = as.numeric(pf[k, ]),
      R2_percent = as.numeric(pf[k, ]) * 100, share_nonnegative = as.numeric(p),
      entropy_factor = factors$entropy_factor[k],
      effective_views_factor = factors$effective_views_factor[k])
  }))
  list(factors = factors, factor_view = fv,
    totals = data.frame(view = VIEWS, R2_total = as.numeric(total),
      R2_total_percent = as.numeric(total) * 100, n_features = as.integer(m$feature_n[VIEWS])),
    min_R2 = min(total), max_R2 = max(total), mean_R2 = mean(total),
    minmax = if (max(positive) > 0) min(positive) / max(positive) else 0,
    entropy = unname(vh["entropy"]), positive_sum = sum(positive),
    negative_total_count = sum(total < 0), negative_factor_view_count = sum(pf < 0))
}

compute_metrics <- function(m, family, scenario, banks) {
  z <- m$Z; g <- m$metadata$grupo
  global <- global_factor_metrics(z, g, banks)
  pw <- pairwise_metrics(z, g, banks)
  pw$q_within_scenario <- stats::p.adjust(pw$p, "BH")
  pw$permanova_q_within_scenario <- stats::p.adjust(pw$permanova_p, "BH")
  si <- silhouette_metrics(z, g, rownames(z))
  di <- dispersion_metrics(z, g, banks)
  disp_rms <- di$table$within_RMS_bias_adjusted
  
  dispersion_balance <- if (
    length(disp_rms) == length(GROUPS) &&
    all(is.finite(disp_rms)) &&
    max(disp_rms) > NUMERICAL_TOL
  ) {
    min(disp_rms) / max(disp_rms)
  } else {
    0
  }
  
  vw <- view_metrics(m)
  signal <- entropy_stats(global$ssb, ncol(z))
  met <- list(family = family, scenario = scenario, n_samples = nrow(z),
    n_NP = sum(g == "NP"), n_FD = sum(g == "FD"), n_SO = sum(g == "SO"),
    n_factors = ncol(z), n_constant_factors = sum(m$zero_factor),
    R2_grupo_global = global$R2, p_grupo_global = global$p,
    silhouette_global = si$global, silhouette_NP = si$by_group["NP"],
    silhouette_FD = si$by_group["FD"], silhouette_SO = si$by_group["SO"],
    silhouette_macro = si$macro, silhouette_worst = si$worst,
    dispersion_balance = dispersion_balance,
    worst_silhouette_group = names(which.min(si$by_group)),
    pairwise_separation_macro = mean(pw$separation), pairwise_separation_worst = min(pw$separation),
    worst_contrast = pw$contrast[which.min(pw$separation)],
    mean_factor_r2_group = mean(global$factors$r2_group),
    median_factor_r2_group = stats::median(global$factors$r2_group),
    max_factor_r2_group = max(global$factors$r2_group), min_factor_q_group = min(global$factors$q_grupo),
    n_factor_group_q010 = sum(global$factors$q_grupo <= 0.10),
    n_factor_group_q005 = sum(global$factors$q_grupo <= 0.05),
    max_group_signal_factor_share = signal["max_share"],
    effective_group_signal_factors = signal["effective"],
    min_view_R2 = vw$min_R2, max_view_R2 = vw$max_R2, mean_view_R2 = vw$mean_R2,
    view_balance_entropy = vw$entropy, view_balance_minmax = vw$minmax,
    mean_factor_view_entropy = mean(vw$factors$entropy_factor),
    median_factor_view_entropy = stats::median(vw$factors$entropy_factor),
    min_factor_view_entropy = min(vw$factors$entropy_factor),
    mean_effective_views = mean(vw$factors$effective_views_factor),
    n_zero_contribution_factors = sum(vw$factors$zero_positive_contribution),
    negative_total_R2_count = vw$negative_total_count,
    negative_factor_view_R2_count = vw$negative_factor_view_count,
    p_dispersion_global = di$p, dispersion_method = di$method,
    baseline = grepl("_qtx00_qpr00_qme00_qcl00$", scenario),
    model_path = m$source, source_md5 = m$source_md5, extraction_backend = m$backend,
    data_storage = m$data_storage, z_geometry = Z_GEOMETRY)
  for (i in seq_len(nrow(pw))) {
    nm <- pw$contrast[i]
    met[[paste0("separation_", nm)]] <- pw$separation[i]
    met[[paste0("p_", nm)]] <- pw$p[i]
    met[[paste0("permanova_R2_", nm)]] <- pw$permanova_R2[i]
    met[[paste0("permanova_p_", nm)]] <- pw$permanova_p[i]
  }
  for (v in VIEWS) {
    met[[paste0("R2_total_", v)]] <- m$view_total[v]
    met[[paste0("n_features_", v)]] <- m$feature_n[v]
  }
  required_metrics <- c(QUALITY_COMPONENTS, paste0("n_features_", VIEWS))
  required_values <- unlist(met[required_metrics], use.names = FALSE)
  met$feasible <- all(is.finite(required_values)) &&
    all(as.numeric(unlist(met[paste0("n_features_", VIEWS)])) > 0) &&
    vw$positive_sum > 0
  met$feasibility_reason <- if (met$feasible) {
    "technical_checks_passed"
  } else {
    "nonfinite_quality_metric_nonpositive_feature_count_or_no_positive_view_R2"
  }
  tag <- function(x) data.frame(family = family, scenario = scenario, x, row.names = NULL)
  list(metrics = as.data.frame(lapply(met, unname), stringsAsFactors = FALSE),
    factor_group = tag(global$factors), pairwise = tag(pw), factor_view = tag(vw$factor_view),
    factor_integration = tag(vw$factors), view_totals = tag(vw$totals),
    silhouette_samples = tag(si$samples), dispersion = tag(di$table),
    model = m, dispersion_note = di$note)
}

rank_scenarios <- function(metrics) {
  valid <- metrics[which(metrics$feasible), , drop = FALSE]
  if (!nrow(valid)) return(valid)

  # q00 se conserva como baseline diagnostico, pero NO participa
  # en la deteccion de meseta ni puede ser seleccionado.
  valid$selection_eligible <- !valid$baseline
  eligible_idx <- which(valid$selection_eligible)
  assert(length(eligible_idx) > 0L,
         "No hay escenarios filtrados elegibles; solo baseline o ninguno")

  # 1) Calidad geometrica/multivista del TRAIN.
  valid$quality_score <- quality_score(valid)
  valid$quality_rank <- NA_integer_
  valid$quality_rank[eligible_idx] <- rank(
    -valid$quality_score[eligible_idx], ties.method = "min"
  )

  plateau_local <- detect_quality_plateau(valid$quality_score[eligible_idx])
  valid$quality_plateau <- FALSE
  valid$quality_plateau[eligible_idx] <- plateau_local$in_plateau
  valid$quality_plateau_cutoff <- plateau_local$cutoff
  valid$quality_plateau_next_score <- if (!is.null(plateau_local$next_score)) {
    plateau_local$next_score
  } else {
    NA_real_
  }
  valid$quality_plateau_largest_drop <- plateau_local$largest_drop
  valid$quality_plateau_split_after <- plateau_local$split_after

  # 2) Cobertura biologica equilibrada entre vistas.
  # Preferentemente relativa al baseline q00; si falta, usa maximo familiar.
  bio <- biological_coverage(valid)
  valid$biological_coverage_score <- bio$score
  valid$biological_coverage_rank <- NA_integer_
  valid$biological_coverage_rank[eligible_idx] <- rank(
    -valid$biological_coverage_score[eligible_idx], ties.method = "min"
  )
  valid$minimum_view_retention <- bio$minimum
  valid$biological_reference_source <- bio$reference_source

  for (v in VIEWS) {
    valid[[paste0("retention_", v)]] <- bio$ratios[, v]
    valid[[paste0("reference_features_", v)]] <- unname(bio$reference[v])
  }

  # 3) Identificar de forma unica el maximo de calidad entre elegibles.
  qopt_local <- eligible_idx[which.max(valid$quality_score[eligible_idx])]
  valid$metric_quality_optimum <- FALSE
  valid$metric_quality_optimum[qopt_local] <- TRUE

  # 4) Seleccion: solo escenarios elegibles; primero meseta de alta calidad;
  # dentro de ella, maxima cobertura biologica; despues proteccion de la peor
  # vista; quality_score actua como desempate posterior.
  ord <- order(
    !valid$selection_eligible,
    !valid$quality_plateau,
    -valid$biological_coverage_score,
    -valid$minimum_view_retention,
    -valid$quality_score,
    -valid$pairwise_separation_worst,
    valid$scenario,
    method = "radix"
  )

  valid <- valid[ord, , drop = FALSE]
  rownames(valid) <- NULL
  valid$final_rank <- seq_len(nrow(valid))

  valid$selected_scenario <- FALSE
  first_eligible <- which(valid$selection_eligible)[1L]
  valid$selected_scenario[first_eligible] <- TRUE

  valid$selection_status <- ifelse(
    !valid$selection_eligible,
    "BASELINE_DIAGNOSTIC_ONLY",
    ifelse(
      valid$quality_plateau,
      "HIGH_QUALITY_PLATEAU",
      "OUTSIDE_HIGH_QUALITY_PLATEAU"
    )
  )
  valid$n_high_quality_plateau <- sum(valid$quality_plateau & valid$selection_eligible)

  valid
}

baseline_comparison <- function(ranking) {
  if (!nrow(ranking) || !any(ranking$baseline)) return(data.frame())

  selected <- ranking[ranking$selected_scenario, , drop = FALSE]
  baseline <- ranking[ranking$baseline, , drop = FALSE]
  keys <- c(
    QUALITY_COMPONENTS,
    "quality_score",
    "biological_coverage_score",
    "minimum_view_retention",
    "pairwise_separation_macro"
  )

  selected_values <- as.numeric(unlist(
    selected[1, keys, drop = FALSE],
    use.names = FALSE
  ))
  baseline_values <- as.numeric(unlist(
    baseline[1, keys, drop = FALSE],
    use.names = FALSE
  ))

  data.frame(
    family = selected$family,
    selected = selected$scenario,
    baseline = baseline$scenario[1],
    metric = keys,
    selected_value = selected_values,
    baseline_value = baseline_values,
    delta_selected_minus_baseline = selected_values - baseline_values,
    stringsAsFactors = FALSE
  )
}

short_scenario <- function(s) {
  out <- sub("^.*_qtx", "tx", s)
  out <- gsub("_qpr", " / pr", out, fixed = TRUE)
  out <- gsub("_qme", " / me", out, fixed = TRUE)
  gsub("_qcl", " / cl", out, fixed = TRUE)
}
plot_style <- function(left = 11, bottom = 4.5) {
  graphics::par(mar = c(bottom, left, 4, 2), bg = "white", las = 1,
                family = "sans", col.axis = "#334155", col.lab = "#172B4D", cex = 0.85)
}
blank_plot <- function(message) {
  graphics::plot.new(); graphics::text(0.5, 0.5, message, cex = 1.1)
}
label_points <- function(xy, labels, colors, cex = 0.68) {
  # Posiciones deterministas; las etiquetas se mueven, los puntos nunca.
  boxes <- list(); h <- graphics::strheight("Mg", cex = cex) * 1.35
  lim <- graphics::par("usr")
  for (i in seq_len(nrow(xy))) {
    w <- graphics::strwidth(labels[i], cex = cex) * 1.1
    candidates <- do.call(rbind, lapply(c(1, 1.6, 2.3, 3.2, 4.5), function(mult) {
      dx <- (w / 2 + h * 0.4) * mult; dy <- h * 1.4 * mult
      cbind(c(0, 0, dx, -dx, dx, -dx, dx, -dx),
            c(dy, -dy, 0, 0, dy, dy, -dy, -dy))
    }))
    chosen <- NULL
    for (j in seq_len(nrow(candidates))) {
      p <- xy[i, ] + candidates[j, ]; box <- c(p[1] - w/2, p[1] + w/2, p[2] - h/2, p[2] + h/2)
      inside <- box[1] >= lim[1] && box[2] <= lim[2] && box[3] >= lim[3] && box[4] <= lim[4]
      collision <- any(vapply(boxes, function(b) box[1] < b[2] && box[2] > b[1] &&
                               box[3] < b[4] && box[4] > b[3], logical(1)))
      over_point <- any(xy[,1] > box[1] - h*.3 & xy[,1] < box[2] + h*.3 &
                          xy[,2] > box[3] - h*.3 & xy[,2] < box[4] + h*.3)
      if (inside && !collision && !over_point) { chosen <- p; break }
    }
    if (is.null(chosen)) {
      chosen <- xy[i, ] + candidates[nrow(candidates), ]
      box <- c(chosen[1] - w/2, chosen[1] + w/2, chosen[2] - h/2, chosen[2] + h/2)
    }
    boxes[[length(boxes) + 1L]] <- box
    if (sum((chosen - xy[i, ])^2) > (2*h)^2) graphics::segments(
      xy[i,1], xy[i,2], chosen[1], chosen[2], col = "#CBD5E1", lwd = 0.6)
    graphics::text(chosen[1], chosen[2], labels = labels[i], cex = cex, col = colors[i])
  }
}
dot_comparison <- function(r, cols, labels, colors, xlab, main, limits = NULL, reference_zero = FALSE) {
  plot_style()
  y <- rev(seq_len(nrow(r))); vals <- as.matrix(r[, cols, drop = FALSE])
  if (is.null(limits)) limits <- range(c(0, vals), finite = TRUE)
  if (diff(limits) == 0) limits <- limits + c(-0.1, 0.1)
  graphics::plot(NA, xlim = limits, ylim = c(0.2, nrow(r) + 1.3), yaxt = "n",
    xlab = xlab, ylab = "", main = main, bty = "n")
  graphics::abline(h = y, col = "#E5E7EB", lwd = 0.6)
  if (reference_zero) graphics::abline(v = 0, col = "#94A3B8", lty = 2)
  graphics::axis(2, at = y, labels = paste0(r$final_rank, ". ", short_scenario(r$scenario)),
                 tick = FALSE, cex.axis = 0.8)
  offsets <- seq(-0.19, 0.19, length.out = length(cols))
  for (j in seq_along(cols)) graphics::points(vals[, j], y + offsets[j], pch = 15 + j,
                                              col = colors[j], cex = 1.0)
  graphics::legend("topright", legend = labels, col = colors, pch = 16:(15 + length(cols)),
                    horiz = TRUE, bty = "n", cex = 0.85)
}

make_plot_functions <- function(ranking, results) {
  r <- utils::head(ranking, TOP_N)
  selected <- ranking[ranking$selected_scenario, , drop = FALSE]
  quality_opt <- ranking[ranking$metric_quality_optimum, , drop = FALSE]
  selected_model <- results[[selected$scenario]]$model
  pal <- c(NP = "#2367A4", FD = "#D56B19", SO = "#7451A8")
  quality_labels <- c(
    "Peor contraste",
    "Peor silhouette",
    "Silhouette macro",
    "Balance dispersion",
    "Balance vistas",
    "Quality score",
    "Cobertura biologica"
  )

  list(
    `01_RANKING` = function() {
      plot_style(left = 5, bottom = 5)

      x <- ranking$quality_score
      y <- ranking$biological_coverage_score
      plateau <- ranking$quality_plateau
      baseline <- ranking$baseline

      xr <- range(x, finite = TRUE)
      yr <- range(y, finite = TRUE)
      xpad <- max(diff(xr) * 0.10, 0.01)
      ypad <- max(diff(yr) * 0.10, 0.01)

      graphics::plot(
        x, y,
        pch = ifelse(baseline, 17, 16),
        col = ifelse(plateau, "#0B806D", "#94A3B8"),
        cex = 1.05,
        xlim = xr + c(-xpad, xpad),
        ylim = pmax(0, yr + c(-ypad, ypad)),
        xlab = "Quality score geometrico/multivista (0-1)",
        ylab = "Cobertura biologica equilibrada (0-1)",
        main = paste(selected$family, "| Calidad TRAIN vs cobertura biologica"),
        bty = "l"
      )

      cutoff <- selected$quality_plateau_cutoff
      graphics::abline(v = cutoff, lty = 3, col = "#0B806D", lwd = 1.2)

      # Optimo puramente geometrico.
      graphics::points(
        quality_opt$quality_score,
        quality_opt$biological_coverage_score,
        pch = 1, cex = 2.0, lwd = 1.8, col = "#172B4D"
      )

      # Escenario finalmente seleccionado.
      graphics::points(
        selected$quality_score,
        selected$biological_coverage_score,
        pch = 8, cex = 2.0, lwd = 2, col = "#B45309"
      )

      ids <- unique(c(
        which(ranking$selected_scenario),
        which(ranking$metric_quality_optimum),
        which(ranking$baseline)
      ))

      # Etiquetar solo las referencias clave para evitar saturar la figura.
      graphics::text(
        ranking$quality_score[ids],
        ranking$biological_coverage_score[ids],
        labels = short_scenario(ranking$scenario[ids]),
        pos = 3, cex = 0.63,
        col = ifelse(ranking$quality_plateau[ids], "#0B806D", "#64748B")
      )

      graphics::legend(
        "bottomleft",
        legend = c(
          "Meseta alta calidad",
          "Fuera de meseta",
          "Seleccionado",
          "Maximo quality score",
          "Baseline q00",
          "Corte empirico de meseta"
        ),
        col = c("#0B806D", "#94A3B8", "#B45309", "#172B4D", "#334155", "#0B806D"),
        pch = c(16, 16, 8, 1, 17, NA),
        lty = c(NA, NA, NA, NA, NA, 3),
        lwd = c(NA, NA, NA, NA, NA, 1.2),
        bty = "n", cex = 0.78
      )

      graphics::mtext(
        sprintf(
          "Seleccion: %s | quality=%.3f | cobertura=%.3f | meseta=%d escenarios",
          short_scenario(selected$scenario), selected$quality_score,
          selected$biological_coverage_score, selected$n_high_quality_plateau
        ),
        side = 1, line = 3.6, cex = 0.72
      )
    },

    `02_HEATMAP_METRICAS` = function() {
      plot_style(bottom = 9)
      heat_keys <- c(QUALITY_COMPONENTS, "quality_score", "biological_coverage_score")
      raw <- as.matrix(r[, heat_keys, drop = FALSE])

      color_values <- sapply(heat_keys, function(k) {
        all_values <- ranking[[k]]
        if (length(unique(all_values)) == 1L) return(rep(0.5, nrow(ranking)))
        (rank(all_values, ties.method = "average") - 1) / (length(all_values) - 1)
      })
      color_values <- matrix(
        color_values,
        nrow = nrow(ranking), ncol = length(heat_keys)
      )[seq_len(nrow(r)), , drop = FALSE]

      graphics::plot(
        NA,
        xlim = c(0.5, length(heat_keys) + 0.5),
        ylim = c(0.5, nrow(r) + 0.5),
        axes = FALSE, xlab = "", ylab = "",
        main = paste(selected$family, "| Calidad y cobertura biologica")
      )

      colors <- grDevices::colorRampPalette(c("#F1F5F9", "#A8D8CF", "#17665B"))(101)
      for (i in seq_len(nrow(r))) for (j in seq_along(heat_keys)) {
        yy <- nrow(r) + 1 - i
        level <- color_values[i, j]
        graphics::rect(
          j - 0.5, yy - 0.5, j + 0.5, yy + 0.5,
          col = colors[1 + round(level * 100)], border = "white"
        )
        graphics::text(
          j, yy, sprintf("%.3f", raw[i, j]), cex = 0.72,
          col = if (level > 0.7) "white" else "#183B3B"
        )
      }

      graphics::axis(
        2, at = rev(seq_len(nrow(r))),
        labels = paste0(
          r$final_rank, ". ", short_scenario(r$scenario),
          ifelse(r$selected_scenario, "  [SEL]", ""),
          ifelse(r$metric_quality_optimum, "  [QMAX]", "")
        ),
        tick = FALSE, cex.axis = 0.72
      )
      graphics::axis(
        1, at = seq_along(heat_keys), labels = quality_labels,
        las = 2, tick = FALSE, cex.axis = 0.78
      )
      graphics::mtext(
        "Cifras = valores reales. Color = posicion relativa dentro de cada metrica y familia.",
        side = 1, line = 7.7, cex = 0.72
      )
    },

    `03_SEPARACION_BALANCE` = function() {
      plot_style(left = 5)
      graphics::plot(
        ranking$view_balance_entropy,
        ranking$pairwise_separation_worst,
        col = ifelse(ranking$quality_plateau, "#0B806D", "#94A3B8"),
        pch = ifelse(ranking$baseline, 17, 16), cex = 1.1,
        xlim = c(0, 1.04),
        ylim = c(0, min(1, max(ranking$pairwise_separation_worst) * 1.2 + 0.03)),
        xlab = "Entropia normalizada entre vistas",
        ylab = "Peor contraste: separacion equilibrada",
        main = paste(selected$family, "| Separacion y balance de vistas"),
        bty = "l"
      )
      ids <- unique(c(
        which(ranking$baseline),
        which(ranking$selected_scenario),
        which(ranking$metric_quality_optimum)
      ))
      graphics::text(
        ranking$view_balance_entropy[ids],
        ranking$pairwise_separation_worst[ids],
        labels = paste0(
          ifelse(ranking$selected_scenario[ids], "SEL ", ""),
          ifelse(ranking$metric_quality_optimum[ids], "QMAX ", ""),
          short_scenario(ranking$scenario[ids])
        ),
        pos = 3, cex = 0.65
      )
      graphics::legend(
        "topleft",
        c("Meseta alta calidad", "Fuera de meseta", "Baseline q00"),
        col = c("#0B806D", "#94A3B8", "#334155"),
        pch = c(16, 16, 17), bty = "n"
      )
    },

    `04_SILHOUETTE_GRUPOS` = function() dot_comparison(
      r, paste0("silhouette_", GROUPS), GROUPS, pal[GROUPS],
      "Silhouette media por grupo",
      paste(selected$family, "| NP, FD y SO tienen igual peso en la macro"),
      c(-1, 1), TRUE
    ),

    `05_SEPARACION_CONTRASTES` = function() dot_comparison(
      r,
      c("separation_NP_FD", "separation_NP_SO", "separation_FD_SO"),
      c("NP-FD", "NP-SO", "FD-SO"),
      c("#D56B19", "#2367A4", "#7451A8"),
      "Separacion equilibrada (0-1)",
      paste(selected$family, "| Los tres contrastes"), c(0, 1)
    ),

    `06_VARIANZA_VISTAS` = function() dot_comparison(
      utils::head(r, 10),
      paste0("R2_total_", VIEWS), c("Tx", "Pr", "Me", "Clinical"),
      c("#2367A4", "#D56B19", "#0B806D", "#7451A8"),
      "Varianza explicada total por vista (proporcion)",
      paste(selected$family, "| Reconstruccion de las cuatro vistas"),
      range(c(0, 1, as.matrix(r[, paste0("R2_total_", VIEWS), drop = FALSE])))
    ),

    `07_FACTOR_VISTA` = function() {
      top <- utils::head(ranking, TOP_FACTOR_PLOTS)
      rows <- bind_rows(lapply(seq_len(nrow(top)), function(i) {
        x <- results[[top$scenario[i]]]$model$factor_view
        data.frame(
          label = paste0("#", top$final_rank[i], " ", short_scenario(top$scenario[i]), " | ", rownames(x)),
          x, check.names = FALSE, row.names = NULL
        )
      }))
      plot_style(left = 15, bottom = 5)
      mat <- as.matrix(rows[, VIEWS, drop = FALSE])
      n <- nrow(mat)
      graphics::plot(
        NA, xlim = c(0.5, 4.5), ylim = c(0.5, n + 0.5), axes = FALSE,
        xlab = "", ylab = "", main = paste(selected$family, "| Factor x vista: R2 en %")
      )
      colors <- grDevices::colorRampPalette(c("#F1F5F9", "#8ECACB", "#155B70"))(101)
      max_color <- max(pmax(mat, 0), NUMERICAL_TOL)
      for (i in seq_len(n)) for (j in 1:4) {
        value <- max(0, mat[i, j]) / max_color
        graphics::rect(
          j - 0.5, n + 1 - i - 0.5, j + 0.5, n + 1 - i + 0.5,
          col = colors[1 + round(value * 100)], border = "white"
        )
        graphics::text(
          j, n + 1 - i, sprintf("%.2f", mat[i, j] * 100), cex = 0.75,
          col = if (value > 0.65) "white" else "#172B4D"
        )
      }
      graphics::axis(2, at = rev(seq_len(n)), labels = rows$label, tick = FALSE, cex.axis = 0.67)
      graphics::axis(1, at = 1:4, labels = VIEWS, tick = FALSE, cex.axis = 0.9)
      graphics::mtext(
        "Escala de color comun a todos los factores mostrados; cifras = R2 porcentual.",
        side = 1, line = 3.5, cex = 0.75
      )
    },

    `08_PROYECCION_Z_GANADOR` = function() {
      plot_style(left = 5, bottom = 6)
      z <- selected_model$Z
      pc <- stats::prcomp(z, center = TRUE, scale. = FALSE)
      explained <- pc$sdev^2 / sum(pc$sdev^2)
      xy <- cbind(pc$x[, 1], if (ncol(pc$x) >= 2L) pc$x[, 2] else 0)
      g <- selected_model$metadata$grupo
      expand <- function(v) {
        rr <- range(v)
        pad <- max(diff(rr) * 0.18, 0.4)
        rr + c(-pad, pad)
      }
      graphics::plot(
        xy, pch = c(NP = 16, FD = 17, SO = 15)[g], col = pal[g], cex = 1.3,
        xlim = expand(xy[, 1]), ylim = expand(xy[, 2]),
        xlab = sprintf("PC1 de Z (%.1f%%)", 100 * explained[1]),
        ylab = if (length(explained) >= 2) sprintf("PC2 de Z (%.1f%%)", 100 * explained[2]) else "K=1: segunda coordenada cero",
        main = paste(selected$family, "| Escenario seleccionado - TRAIN"),
        bty = "l", asp = 1
      )
      label_points(xy, selected_model$metadata$sample, pal[g])
      graphics::legend(
        "topright", legend = c("NP (6)", "FD (3)", "SO (7)"),
        col = pal[GROUPS], pch = c(16, 17, 15), bty = "n"
      )
      graphics::mtext(short_scenario(selected$scenario), side = 3, line = 0.3, cex = 0.8)
      graphics::mtext(
        "Proyeccion descriptiva; la seleccion combina meseta de calidad TRAIN y cobertura biologica.",
        side = 1, line = 4.8, cex = 0.75
      )
    }
  )
}

write_plots <- function(ranking, results, out) {
  if (!MAKE_PLOTS || !nrow(ranking)) return(data.frame())
  pdir <- file.path(out, "plots"); dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  functions <- make_plot_functions(ranking, results)
  audit <- list()
  for (nm in names(functions)) for (format in c("pdf", "png")) {
    path <- file.path(pdir, paste0(nm, ".", format))
    nr <- if (nm == "07_FACTOR_VISTA") sum(vapply(utils::head(ranking$scenario, TOP_FACTOR_PLOTS),
      function(s) ncol(results[[s]]$model$Z), integer(1))) else min(TOP_N, nrow(ranking))
    height <- if (nm %in% c("01_RANKING", "03_SEPARACION_BALANCE", "08_PROYECCION_Z_GANADOR")) 7.5 else max(6, nr * 0.34 + 3)
    device <- NULL
    previous_device <- grDevices::dev.cur()
    warnings <- character()
    err <- tryCatch(withCallingHandlers({
      if (format == "pdf") grDevices::pdf(path, width = 12, height = height, useDingbats = FALSE) else {
        grDevices::png(path, width = 12, height = height, units = "in", res = 140,
                       type = if (capabilities("cairo")) "cairo" else getOption("bitmapType", "Xlib"))
      }
      assert(grDevices::dev.cur() > 1L && grDevices::dev.cur() != previous_device,
             paste("No se pudo abrir el dispositivo", format, "; se evita crear Rplots.pdf fuera de downstream"))
      device <- grDevices::dev.cur()
      functions[[nm]]()
      grDevices::dev.off(device); device <- NULL
      assert(file.exists(path) && file.info(path)$size > 0, "El dispositivo no genero un archivo valido")
      ""
    }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) conditionMessage(e))
    if (!is.null(device)) try(grDevices::dev.off(device), silent = TRUE)
    audit[[length(audit) + 1L]] <- data.frame(plot = nm, format = format,
      status = if (nzchar(err)) "ERROR" else if (length(warnings)) "OK_WITH_WARNING" else "OK",
      error = err, warnings = paste(unique(warnings), collapse = " | "), path = path)
  }
  bind_rows(audit)
}

summary_text <- function(family, ranking, audit) {
  header <- c(
    paste("METATEST |", family, "| SELECCION EXCLUSIVA TRAIN"),
    paste("Directorios auditados:", nrow(audit)),
    paste("Escenarios factibles:", nrow(ranking)),
    ""
  )

  if (!nrow(ranking)) {
    return(c(header, "SIN SELECCION: no hay escenarios factibles. Revisar SCENARIO_AUDIT.csv."))
  }

  selected <- ranking[ranking$selected_scenario, , drop = FALSE]
  qopt <- ranking[ranking$metric_quality_optimum, , drop = FALSE]
  plateau <- ranking[ranking$quality_plateau, , drop = FALSE]

  text <- c(
    header,
    paste("ESCENARIO SELECCIONADO:", selected$scenario),
    paste("OPTIMO PURAMENTE GEOMETRICO:", qopt$scenario),
    paste("Estado:", selected$selection_status),
    paste("Modelo original:", selected$model_path),
    paste("MD5:", selected$source_md5),
    paste("TRAIN:", selected$n_samples, "muestras; NP=6, FD=3, SO=7."),
    paste("Factores:", selected$n_factors, "; geometria:", Z_GEOMETRY, "; se usan todos."),
    "",
    paste("Quality score seleccionado / maximo:", fmt(selected$quality_score), "/", fmt(qopt$quality_score)),
    paste("Cobertura biologica seleccionada:", fmt(selected$biological_coverage_score)),
    paste("Peor retencion entre vistas:", fmt(selected$minimum_view_retention)),
    paste("Meseta empirica de alta calidad:", nrow(plateau), "escenarios."),
    paste("Corte quality_score:", fmt(selected$quality_plateau_cutoff)),
    paste("Mayor salto que define la meseta:", fmt(selected$quality_plateau_largest_drop)),
    "",
    "COMPONENTES DE CALIDAD DEL SELECCIONADO:",
    paste("Peor contraste:", fmt(selected$pairwise_separation_worst)),
    paste("Peor silhouette:", fmt(selected$silhouette_worst)),
    paste("Silhouette macro:", fmt(selected$silhouette_macro)),
    paste("Balance de dispersion:", fmt(selected$dispersion_balance)),
    paste("Balance de vistas (entropia):", fmt(selected$view_balance_entropy)),
    "",
    "RETENCION BIOLOGICA POR VISTA:",
    paste("Tx / Pr / Me / Clinical:", paste(
      fmt(as.numeric(unlist(
        selected[1, paste0("retention_", VIEWS), drop = FALSE],
        use.names = FALSE
      ))),
      collapse = " / "
    )),
    paste("Referencia de cobertura:", selected$biological_reference_source),
    "",
    paste("Separacion NP-FD / NP-SO / FD-SO:", paste(fmt(c(
      selected$separation_NP_FD,
      selected$separation_NP_SO,
      selected$separation_FD_SO
    )), collapse = " / ")),
    paste("Silhouette NP / FD / SO:", paste(fmt(c(
      selected$silhouette_NP,
      selected$silhouette_FD,
      selected$silhouette_SO
    )), collapse = " / ")),
    paste("R2 total Tx / Pr / Me / Clinical (%):", paste(
      fmt(100 * as.numeric(unlist(
        selected[1, paste0("R2_total_", VIEWS), drop = FALSE],
        use.names = FALSE
      ))),
      collapse = " / "
    )),
    "",
    "REGLA DE SELECCION:",
    "1) Calcular quality_score como media geometrica de separacion peor, silhouette peor, silhouette macro, balance de dispersion y balance de vistas.",
    "2) Detectar empiricamente la meseta de alta calidad mediante el mayor salto del quality_score ordenado.",
    "3) Dentro de esa meseta, seleccionar la maxima cobertura biologica multivista (media geometrica de retenciones).",
    "4) Desempatar por peor retencion de vista, quality_score, peor separacion y nombre de escenario.",
    "",
    "LIMITACIONES:",
    "La meseta es empirica; no es un intervalo de confianza ni prueba de equivalencia.",
    "Con N=16 y FD=3, la seleccion es exploratoria y condicional a los modelos MOFA finales dados.",
    "La cobertura biologica mide retencion de features por vista; no demuestra mayor varianza biologica funcional.",
    "Los p/q son diagnosticos in-sample; no validan superioridad poblacional.",
    "No se estima incertidumbre del ajuste MOFA ni estabilidad entre semillas.",
    "No se ha leido, proyectado ni evaluado ningun archivo TEST."
  )

  if (any(ranking$baseline)) {
    base <- ranking[ranking$baseline, , drop = FALSE]
    text <- c(
      text, "",
      paste("Baseline q00:", base$scenario[1]),
      paste("Quality score baseline:", fmt(base$quality_score[1])),
      paste("Cobertura baseline:", fmt(base$biological_coverage_score[1]))
    )
  }

  if (selected$negative_total_R2_count + selected$negative_factor_view_R2_count > 0) {
    text <- c(text, "Hay R2 negativos: se conservan en tablas; solo para participaciones/entropias se truncan en cero.")
  }
  if (selected$silhouette_FD < 0 || selected$silhouette_worst < 0) {
    text <- c(text, "ADVERTENCIA: al menos un grupo tiene silhouette negativa.")
  }

  text
}

run_family <- function(family, input, output, manifest, banks, parameters) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  audit <- list(); results <- list()
  directories <- if (dir.exists(input)) sort(list.dirs(input, full.names = TRUE, recursive = FALSE)) else character()
  if (!length(directories)) audit[[1L]] <- data.frame(family = family, scenario = NA_character_,
    status = "NO_SCENARIO_DIRECTORIES", error = paste("Carpeta inexistente o vacia:", input), model_path = NA_character_)
  for (directory in directories) {
    scenario <- basename(directory)
    path <- file.path(directory, "integracion", "rds", paste0("modelo_mofa_final_", scenario, ".rds"))
    message(sprintf("[%s] %s", family, scenario))
    started <- proc.time()["elapsed"]; warnings <- character()
    result <- tryCatch(withCallingHandlers({
      assert(grepl(FAMILY_PATTERNS[[family]], scenario), "Nombre de escenario incompatible con familia; no se mezcla")
      assert(file.exists(path), "No existe el modelo final en la ruta exacta esperada")
      m <- extract_model(path, manifest)
      compute_metrics(m, family, scenario, banks)
    }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
    failed <- inherits(result, "error")
    aud <- data.frame(family = family, scenario = scenario, model_path = path,
      status = if (failed) "REJECTED" else if (!result$metrics$feasible) "INFEASIBLE" else "OK",
      error = if (failed) conditionMessage(result) else "",
      warnings = paste(unique(warnings), collapse = " | "),
      seconds = unname(proc.time()["elapsed"] - started), stringsAsFactors = FALSE)
    if (!failed) {
      results[[scenario]] <- result
      aud$source_md5 <- result$model$source_md5
      aud$z_reordered_to_manifest <- result$model$reordered
      aud$backend <- result$model$backend; aud$data_storage <- result$model$data_storage
      aud$train_id_check <- "EXACT_MANIFEST_MATCH_METADATA_Z_AND_ALL_VIEWS"
      aud$dispersion_note <- result$dispersion_note
    } else message("  RECHAZADO: ", conditionMessage(result))
    audit[[length(audit) + 1L]] <- aud
    write_csv(bind_rows(audit), file.path(output, "SCENARIO_AUDIT.csv"))
  }
  aud <- bind_rows(audit)
  raw_metrics <- bind_rows(lapply(results, `[[`, "metrics"))
  ranking <- if (nrow(raw_metrics)) {
    tryCatch(
      rank_scenarios(raw_metrics),
      error = function(e) {
        stop(
          paste0("Fallo en rank_scenarios(): ", conditionMessage(e)),
          call. = FALSE
        )
      }
    )
  } else {
    data.frame()
  }
  all_metrics <- raw_metrics
  if (nrow(ranking)) {
    add <- setdiff(names(ranking), names(raw_metrics))
    idx <- match(all_metrics$scenario, ranking$scenario)
    for (nm in add) all_metrics[[nm]] <- ranking[[nm]][idx]
  }
  rejected <- aud[which(aud$status == "REJECTED"), , drop = FALSE]
  if (nrow(rejected)) all_metrics <- bind_rows(list(all_metrics, data.frame(
    family = family, scenario = rejected$scenario, model_path = rejected$model_path,
    feasible = FALSE, feasibility_reason = rejected$error)))
  tables <- list(
    ALL_SCENARIOS_METRICS = all_metrics,
    FACTOR_GROUP_ASSOCIATION = bind_rows(lapply(results, `[[`, "factor_group")),
    FACTOR_VIEW_VARIANCE = bind_rows(lapply(results, `[[`, "factor_view")),
    FACTOR_VIEW_INTEGRATION = bind_rows(lapply(results, `[[`, "factor_integration")),
    VIEW_TOTAL_VARIANCE = bind_rows(lapply(results, `[[`, "view_totals")),
    PAIRWISE_GROUP_SEPARATION = bind_rows(lapply(results, `[[`, "pairwise")),
    SAMPLE_SILHOUETTE = bind_rows(lapply(results, `[[`, "silhouette_samples")),
    GROUP_DISPERSION = bind_rows(lapply(results, `[[`, "dispersion")),
    SCENARIO_AUDIT = aud,
    SELECTION_RANKING = ranking,
    TOP3_SELECTED_RANKING = utils::head(ranking, 3),
    HIGH_QUALITY_PLATEAU = if (nrow(ranking)) ranking[ranking$quality_plateau, , drop = FALSE] else ranking,
    SELECTED_SCENARIO = if (nrow(ranking)) ranking[ranking$selected_scenario, , drop = FALSE] else ranking,
    METRIC_QUALITY_OPTIMUM = if (nrow(ranking)) ranking[ranking$metric_quality_optimum, , drop = FALSE] else ranking,
    BASELINE_COMPARISON = baseline_comparison(ranking)
  )

  if (nrow(tables$FACTOR_GROUP_ASSOCIATION)) {
    tables$FACTOR_GROUP_ASSOCIATION$q_grupo_family_diagnostic <-
      stats::p.adjust(tables$FACTOR_GROUP_ASSOCIATION$p_grupo, "BH")
  }
  if (nrow(tables$PAIRWISE_GROUP_SEPARATION)) {
    tables$PAIRWISE_GROUP_SEPARATION$q_family_diagnostic <-
      stats::p.adjust(tables$PAIRWISE_GROUP_SEPARATION$p, "BH")
  }

  plots <- tryCatch(
    write_plots(ranking, results, output),
    error = function(e) data.frame(
      plot = "PLOT_SETUP", status = "ERROR", error = conditionMessage(e)
    )
  )
  tables$PLOT_AUDIT <- plots

  for (nm in names(tables)) {
    write_csv(tables[[nm]], file.path(output, paste0(nm, ".csv")))
  }

  writeLines(
    summary_text(family, ranking, aud),
    file.path(output, "SELECTED_SCENARIO_SUMMARY.txt"),
    useBytes = TRUE
  )

  extracted <- lapply(results, function(x) x$model)
  saveRDS(
    list(
      tables = tables,
      parameters = parameters,
      family = family,
      input = input,
      extracted_models = extracted,
      permutation_bank_path = file.path(dirname(output), "PERMUTATION_BANKS.rds")
    ),
    file.path(output, "SCENARIO_SELECTION_OBJECT.rds"),
    version = 3
  )

  plot_error <- nrow(plots) > 0L && any(plots$status != "OK")
  selected <- if (nrow(ranking)) ranking[ranking$selected_scenario, , drop = FALSE] else ranking

  list(
    selected = selected,
    audit = aud,
    status = data.frame(
      family = family,
      input = input,
      n_directories = length(directories),
      n_valid = nrow(ranking),
      n_rejected = sum(aud$status == "REJECTED"),
      n_infeasible = sum(aud$status == "INFEASIBLE"),
      n_high_quality_plateau = if (nrow(ranking)) sum(ranking$quality_plateau) else 0,
      status = if (!nrow(ranking)) {
        "NO_SELECTION"
      } else if (plot_error || any(aud$status != "OK")) {
        "COMPLETE_WITH_ISSUES"
      } else {
        "COMPLETE"
      },
      output = output
    )
  )
}

parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE); i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key %in% c("--help", "-h")) {
      cat("Rscript select_best_upstream_scenario_downstream.R\n",
        "  --upstream DIR      upstream_comun o una de sus carpetas de familia\n",
        "  --downstream DIR    salida nueva dentro de un directorio downstream\n",
        "  --families P1,P2,P3 familias independientes (por defecto las tres)\n",
        "  --train-manifest CSV  lista cerrada sample,grupo\n",
        "  --seed N --n-perm N --top-n N --z-geometry unit_variance|native\n",
        "  --no-plots\n", sep = "")
      return(FALSE)
    }
    if (key == "--no-plots") { MAKE_PLOTS <<- FALSE; i <- i + 1L; next }
    allowed <- c("--upstream", "--downstream", "--families", "--train-manifest",
                 "--seed", "--n-perm", "--top-n", "--z-geometry")
    assert(key %in% allowed && i < length(args), paste("Argumento invalido/incompleto:", key))
    val <- args[i + 1L]
    switch(key, `--upstream` = { UPSTREAM_DIR <<- val }, `--downstream` = { DOWNSTREAM_DIR <<- val },
      `--families` = { FAMILIES <<- strsplit(val, ",", fixed = TRUE)[[1]] },
      `--train-manifest` = { TRAIN_MANIFEST_FILE <<- val }, `--seed` = { SEED <<- as.integer(val) },
      `--n-perm` = { N_PERM <<- as.integer(val) }, `--top-n` = { TOP_N <<- as.integer(val) },
      `--z-geometry` = { Z_GEOMETRY <<- val })
    i <- i + 2L
  }
  TRUE
}
main <- function() {
  if (!parse_cli()) return(invisible(NULL))
  options(stringsAsFactors = FALSE, digits = 17, warn = 1)
  assert(getRversion() >= "4.1", "Se requiere R >= 4.1")
  assert(length(SEED) == 1L && is.finite(SEED) && SEED >= 0, "SEED invalida")
  assert(length(N_PERM) == 1L && is.finite(N_PERM) && N_PERM >= 1L, "N_PERM debe ser entero positivo")
  assert(length(TOP_N) == 1L && is.finite(TOP_N) && TOP_N >= 1, "TOP_N invalido")
  assert(Z_GEOMETRY %in% c("unit_variance", "native"), "Geometria Z no reconocida")
  assert(length(QUALITY_COMPONENTS) == 5L && !anyDuplicated(QUALITY_COMPONENTS),
         "QUALITY_COMPONENTS debe contener cinco metricas unicas")
  assert(!anyDuplicated(FAMILIES) && length(FAMILIES) > 0 && all(FAMILIES %in% names(FAMILY_DIRS)), "Familias invalidas")
  manifest <- if (nzchar(TRAIN_MANIFEST_FILE)) utils::read.csv(path.expand(TRAIN_MANIFEST_FILE),
    stringsAsFactors = FALSE, check.names = FALSE) else TRAIN_MANIFEST
  manifest <- validate_manifest(manifest)
  upstream <- canonical_path(UPSTREAM_DIR)
  assert(dir.exists(upstream), paste("UPSTREAM_DIR no existe:", upstream))
  direct_family <- names(FAMILY_DIRS)[FAMILY_DIRS == basename(upstream)]
  if (length(direct_family)) { FAMILIES <<- direct_family; root <- dirname(upstream) } else root <- upstream
  downstream <- if (nzchar(DOWNSTREAM_DIR)) canonical_path(DOWNSTREAM_DIR) else
    file.path(root, "downstream", "seleccion_mejor_escenario_train")
  assert("downstream" %in% strsplit(downstream, "/", fixed = TRUE)[[1]],
         "DOWNSTREAM_DIR debe estar dentro de una carpeta denominada downstream")
  for (d in FAMILY_DIRS) assert(!is_inside(downstream, canonical_path(file.path(root, d))),
    "La salida no puede estar dentro de una carpeta de escenarios upstream")
  stamp <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_pid", Sys.getpid())
  output <- file.path(downstream, stamp)
  assert(!dir.exists(output), "La carpeta de corrida ya existe; no se sobrescribe")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  assert(dir.exists(output), "No se pudo crear la carpeta downstream")
  packages <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "MOFA2", "vegan", "cluster", "data.table", "ggplot2")
  parameters <- list(script_version = "2.0.1", date_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    R_version = R.version.string, package_versions = setNames(lapply(packages, pkg_version), packages),
    UPSTREAM_DIR = upstream, DOWNSTREAM_DIR = downstream, run_output = output, FAMILIES = FAMILIES,
    SEED = SEED, N_PERM = N_PERM, TOP_N = TOP_N, TOP_FACTOR_PLOTS = TOP_FACTOR_PLOTS,
    Z_GEOMETRY = Z_GEOMETRY, MOFA_R2_DIVISOR = MOFA_R2_DIVISOR,
    QUALITY_COMPONENTS = QUALITY_COMPONENTS, QUALITY_EPS = QUALITY_EPS,
    selection_rule = "empirical_high_quality_plateau_then_max_balanced_biological_coverage",
    plateau_rule = "largest_drop_in_ordered_quality_score",
    biological_coverage = "geometric_mean_of_per_view_feature_retention_relative_to_q00",
    TRAIN_MANIFEST = manifest, TRAIN_MANIFEST_FILE = TRAIN_MANIFEST_FILE,
    train_provenance = "User-declared TRAIN, exact IDs from supplied reference RDS; no split inference from labels",
    primary_effect = "equal_prior_centroid_R2_with_unbiased_within_group_trace",
    permutation_assumption = "independent exchangeable samples under the null",
    uncertainty = "empirical plateau, not statistical equivalence; no MOFA refit uncertainty")
  dput(parameters, file = file.path(output, "PARAMETERS.R"))
  writeLines(capture.output(utils::sessionInfo()), file.path(output, "SESSION_INFO.txt"))
  write_csv(manifest, file.path(output, "TRAIN_MANIFEST_USED.csv"))
  script_arg <- commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))]
  if (length(script_arg)) {
    script_path <- sub("^--file=", "", script_arg[1])
    if (file.exists(script_path)) file.copy(script_path, file.path(output, "SCRIPT_USED.R"), overwrite = FALSE)
  }
  banks <- make_banks(manifest, N_PERM, SEED)
  saveRDS(banks, file.path(output, "PERMUTATION_BANKS.rds"), version = 3)
  res <- list()
  for (family in FAMILIES) {
    out <- file.path(output, family)
    res[[family]] <- tryCatch(run_family(family, file.path(root, FAMILY_DIRS[[family]]),
       out, manifest, banks, parameters), error = function(e) {
      dir.create(out, recursive = TRUE, showWarnings = FALSE)
      writeLines(conditionMessage(e), file.path(out, "FAMILY_ERROR.txt"))
      message("ERROR EN FAMILIA ", family, ": ", conditionMessage(e))
      list(selected = data.frame(), audit = data.frame(family = family, status = "FAMILY_ERROR", error = conditionMessage(e)),
        status = data.frame(family = family, status = "FAMILY_ERROR", output = out))
    })
    write_csv(bind_rows(lapply(res, `[[`, "status")), file.path(output, "RUN_STATUS.csv"))
  }
  selected <- bind_rows(lapply(res, `[[`, "selected"))
  write_csv(selected, file.path(output, "SELECTED_SCENARIOS_BY_ANALYSIS.csv"))
  write_csv(bind_rows(lapply(res, `[[`, "audit")), file.path(output, "SCENARIO_AUDIT_ALL_FAMILIES.csv"))

  writeLines(
    c(
      "Cada fila de SELECTED_SCENARIOS_BY_ANALYSIS.csv procede de una seleccion independiente por familia.",
      "No hay ranking conjunto ni ganador global entre P1, P2 y P3.",
      paste("Salida:", output),
      "Leer RUN_STATUS.csv y cada SELECTED_SCENARIO_SUMMARY.txt.",
      "La seleccion se hace en TRAIN en dos etapas: meseta empirica de alta calidad y, dentro de ella, maxima cobertura biologica equilibrada.",
      "METRIC_QUALITY_OPTIMUM.csv conserva el escenario que maximiza quality_score aunque no sea el seleccionado.",
      "HIGH_QUALITY_PLATEAU.csv contiene los escenarios del regimen empirico de alta calidad.",
      "TEST no se utiliza en ninguna etapa de esta seleccion."
    ),
    file.path(output, "READ_ME_FIRST.txt")
  )

  cat("\nSELECCION DOWNSTREAM TERMINADA\n", output, "\n", sep = "")
  if (nrow(selected)) {
    print(
      selected[, c(
        "family", "scenario", "quality_score",
        "biological_coverage_score", "n_high_quality_plateau"
      ), drop = FALSE],
      row.names = FALSE
    )
  }

  statuses <- bind_rows(lapply(res, `[[`, "status"))$status
  if (any(statuses %in% c("NO_SELECTION", "FAMILY_ERROR"))) {
    stop(
      "Alguna familia solicitada no produjo seleccion. Las otras salidas estan guardadas; revisar RUN_STATUS.csv.",
      call. = FALSE
    )
  }
  if (any(statuses == "COMPLETE_WITH_ISSUES")) message("Hay incidencias: revisar RUN_STATUS.csv, SCENARIO_AUDIT.csv y PLOT_AUDIT.csv.")
  invisible(output)
}
if (sys.nframe() == 0L) main()

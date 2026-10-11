#!/usr/bin/env Rscript

# ==============================================================================
# METATEST — 08_define_biological_programs.R
# ==============================================================================
# PURPOSE
#   Define reproducible transcriptomic biological programs from the already
#   frozen/audited outputs of scripts 01–07, then attach multimodal context.
#
# THIS SCRIPT DOES NOT:
#   - retrain MOFA;
#   - reselect scenarios;
#   - recompute enrichment significance;
#   - use TEST or pooled data to define programs;
#   - select a biomarker shortlist;
#   - require a minimum representation by omic view;
#   - rank P1/P2/P3;
#   - create final figures.
#
# CORE PRINCIPLES
#   1) MOFA factor sign is arbitrary. Phenotypic side is reconstructed from the
#      frozen class<->factor posterior audit for each family independently.
#   2) Primary TX-program definition uses only formal pathway membership among
#      transcriptomic features observed in the frozen MOFA universe.
#   3) Driver genes support interpretation but do not create pathway clusters.
#   4) PR/ME/CL characterize already-defined TX programs; they never determine
#      pathway significance or the initial TX-program partition.
#   5) TEST is transportability audit only. Pooled n=26 is descriptive only.
#   6) P2 vs P3 is treated as sensitivity to clinical specification with/without
#      FLI, not as a single-variable causal ablation of FLI.
#
# PROGRAM DEFINITION
#   - Work within each family x inferred factor role x phenotypic side.
#   - Pairwise redundancy is computed from observed pathway members:
#       Jaccard, overlap coefficient, shared genes, hypergeometric overlap test.
#   - BH-adjusted overlap q <= PROGRAM_OVERLAP_Q defines graph edges.
#   - Weighted Louvain communities define the primary, non-predefined number of
#     programs. Singletons remain valid one-pathway programs.
#   - Cross-community significant edges create SECONDARY affiliations and identify
#     bridge pathways; they do not erase the primary partition.
#   - Representative pathway = structural medoid (highest mean structural
#     similarity to primary program members), with deterministic tie-breaking.
#
# IMPORTANT ABOUT G_annotated vs G_observed
#   Upstream 03 stores the formal MSigDB membership restricted to transcriptomic
#   features actually observed/mapped in the frozen MOFA model. This script uses
#   that frozen G_observed for program definition. It does NOT re-query MSigDB,
#   avoiding version drift. Genes not belonging to any primary pathway are kept
#   explicitly in TX_FEATURE_PRIMARY_PATHWAY_COVERAGE.csv.
# ============================================================================== 

options(stringsAsFactors = FALSE)

# ------------------------------- configuration -------------------------------
UPSTREAM_DIR <- Sys.getenv(
  "UPSTREAM_DIR",
  unset = path.expand("~/Documents/PhD/metatest/results_upstream_local_20260908/comparacion_modelamientos/upstream_comun")
)

SELECTION_RUN <- Sys.getenv("SELECTION_RUN", unset = "")
FAMILIES <- trimws(strsplit(Sys.getenv("BIO_PROGRAM_FAMILIES", unset = "P1,P2,P3"), ",", fixed = TRUE)[[1]])
FAMILIES <- FAMILIES[nzchar(FAMILIES)]

# Reuse the already established 0.10 FDR convention rather than introduce a
# different post-hoc alpha solely for redundancy reduction. It remains explicit
# and can be overridden only deliberately through the environment.
PROGRAM_OVERLAP_Q <- as.numeric(Sys.getenv("BIO_PROGRAM_OVERLAP_Q", unset = "0.10"))
PROGRAM_GRAPH_SEED <- as.integer(Sys.getenv("BIO_PROGRAM_GRAPH_SEED", unset = "20260925"))
NUM_TOL <- 1e-12

# Existing upstream definition: pathway primary at q_both_sides <= 0.10.
PRIMARY_PATHWAY_REQUIRED_FILE <- "04B_MOFA_NATIVE_ENRICHMENT_SIGNIFICANT_TWO_SIDES_Q10.csv"

# ------------------------------- dependencies --------------------------------
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop(
    paste0(
      "Package 'igraph' is required for reproducible graph community detection.\n",
      "Install/reinstall it before running this script. No fallback clustering is used,\n",
      "because changing the clustering algorithm would change the analytical decision."
    ),
    call. = FALSE
  )
}

# --------------------------------- helpers -----------------------------------
assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}

canonical_path <- function(p, must_work = TRUE) {
  normalizePath(path.expand(p), winslash = "/", mustWork = must_work)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x) || !is.data.frame(x) || !ncol(x)) {
    x <- data.frame(status = character(), stringsAsFactors = FALSE)
  }
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

read_csv_checked <- function(path, required = character()) {
  assert(file.exists(path), paste("Missing required file:", path))
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(required, names(x))
  assert(!length(missing), paste("Missing columns in", basename(path), ":", paste(missing, collapse = ", ")))
  x
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes")] <- TRUE
  out[y %in% c("false", "f", "0", "no")] <- FALSE
  out
}

collapse_unique <- function(x, sep = ";") {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return(NA_character_)
  paste(sort(x), collapse = sep)
}

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}

n_true <- function(x) {
  x <- as_logical_safe(x)
  sum(x %in% TRUE, na.rm = TRUE)
}

# --------------------------- locate pipeline runs ----------------------------
find_latest_selection_run <- function(upstream) {
  base <- file.path(upstream, "downstream", "seleccion_mejor_escenario_train")
  assert(dir.exists(base), paste("Selection root does not exist:", base))
  dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "SELECTED_SCENARIOS_BY_ANALYSIS.csv"))]
  assert(length(dirs) > 0L, paste("No complete selection runs found in", base))
  info <- file.info(dirs)
  canonical_path(dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L], TRUE)
}

find_bidirectional_run <- function(selection_run, family) {
  root <- file.path(selection_run, family, "BIDIRECTIONAL_FRAMEWORK")
  assert(dir.exists(root), paste("No BIDIRECTIONAL_FRAMEWORK for", family))
  latest <- file.path(root, "LATEST_RUN.txt")
  if (file.exists(latest)) {
    p <- trimws(readLines(latest, warn = FALSE)[1L])
    if (nzchar(p) && dir.exists(p) && file.exists(file.path(p, "rds", "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"))) {
      return(canonical_path(p, TRUE))
    }
  }
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "rds", "BIDIRECTIONAL_FRAMEWORK_OBJECT.rds"))]
  assert(length(dirs) > 0L, paste("No complete BIDIRECTIONAL_FRAMEWORK run for", family))
  info <- file.info(dirs)
  canonical_path(dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L], TRUE)
}

find_latest_subrun <- function(parent, folder, required_file) {
  root <- file.path(parent, folder)
  assert(dir.exists(root), paste("Missing", folder, "under", parent))
  latest <- file.path(root, "LATEST_RUN.txt")
  if (file.exists(latest)) {
    p <- trimws(readLines(latest, warn = FALSE)[1L])
    if (nzchar(p) && dir.exists(p) && file.exists(file.path(p, required_file))) {
      return(canonical_path(p, TRUE))
    }
  }
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, required_file))]
  assert(length(dirs) > 0L, paste("No complete", folder, "run found under", parent))
  info <- file.info(dirs)
  canonical_path(dirs[order(info$mtime, basename(dirs), decreasing = TRUE)][1L], TRUE)
}

resolve_family_runs <- function(selection_run, family) {
  bidir <- find_bidirectional_run(selection_run, family)
  list(
    family = family,
    bidirectional = bidir,
    latent = find_latest_subrun(bidir, "LATENT_SUPERVISED_AUDIT", file.path("tables", "01_FACTOR_ORIENTATION_AND_SCALING.csv")),
    enrichment = find_latest_subrun(bidir, "MOFA_NATIVE_ENRICHMENT", file.path("tables", PRIMARY_PATHWAY_REQUIRED_FILE)),
    class_pathway = find_latest_subrun(bidir, "CLASS_PATHWAY_BRIDGE", file.path("tables", "01_FACTOR_SIDE_CLASS_POSTERIOR.csv")),
    crossomic = find_latest_subrun(bidir, "CROSSOMIC_BRIDGE", file.path("rds", "CROSSOMIC_BRIDGE_OBJECT.rds")),
    train = find_latest_subrun(bidir, "OBSERVED_TRAIN_CORROBORATION", file.path("rds", "OBSERVED_TRAIN_CORROBORATION_OBJECT.rds")),
    test = find_latest_subrun(bidir, "TEST_TRANSPORTABILITY_AUDIT", file.path("rds", "TEST_TRANSPORTABILITY_AUDIT_OBJECT.rds"))
  )
}

# ----------------------- reproducible phenotype orientation -----------------
infer_factor_roles_and_sides <- function(side_post, family) {
  req <- c(
    "factor", "oriented_side", "p_side_NP", "p_side_FD", "p_side_SO",
    "assigned_class_nominal", "p_assigned_class", "assignment_strong"
  )
  missing <- setdiff(req, names(side_post))
  assert(!length(missing), paste("Side posterior incomplete for", family, ":", paste(missing, collapse = ", ")))

  side_post$assignment_strong <- as_logical_safe(side_post$assignment_strong)
  factors <- unique(as.character(side_post$factor))
  out <- list()
  rr <- 1L

  for (f in factors) {
    z <- side_post[side_post$factor == f, , drop = FALSE]
    assert(nrow(z) == 2L, paste("Expected exactly two sides for", family, f))
    assert(setequal(z$oriented_side, c("positive", "negative")), paste("Unexpected sides for", family, f))

    strong_fd <- which(z$assignment_strong %in% TRUE & z$assigned_class_nominal == "FD")
    strong_np <- which(z$assignment_strong %in% TRUE & z$assigned_class_nominal == "NP")
    strong_so <- which(z$assignment_strong %in% TRUE & z$assigned_class_nominal == "SO")

    role <- "UNRESOLVED_CLASS_AXIS"
    labels <- rep("ambiguous_side", nrow(z))

    # Priority 1: one unambiguous FD side => FD vs non-FD axis.
    if (length(strong_fd) == 1L) {
      role <- "FD_NONFD_AXIS"
      labels[strong_fd] <- "FD_side"
      labels[-strong_fd] <- "nonFD_side"
    } else if (length(strong_np) == 1L && length(strong_so) == 1L && strong_np != strong_so) {
      # Priority 2: opposite strong NP/SO sides => NP vs SO axis.
      role <- "NP_SO_AXIS"
      labels[strong_np] <- "NP_side"
      labels[strong_so] <- "SO_side"
    } else {
      # Generic fallback remains reproducible and does not force a phenotype.
      idx <- which(z$assignment_strong %in% TRUE)
      if (length(idx)) {
        labels[idx] <- paste0(z$assigned_class_nominal[idx], "_side")
      }
    }

    z$family <- family
    z$factor_role <- role
    z$phenotypic_side <- labels
    z$orientation_rule <- if (role == "FD_NONFD_AXIS") {
      "one_strong_FD_side_then_opposite_nonFD"
    } else if (role == "NP_SO_AXIS") {
      "opposite_strong_NP_and_SO_sides"
    } else {
      "generic_strong_assignment_else_ambiguous"
    }

    out[[rr]] <- z
    rr <- rr + 1L
  }

  do.call(rbind, out)
}

# --------------------------- pathway membership -----------------------------
load_primary_pathways <- function(runs, family, side_map) {
  p <- read_csv_checked(
    file.path(runs$enrichment, "tables", PRIMARY_PATHWAY_REQUIRED_FILE),
    required = c(
      "family", "library", "pathway", "factor", "oriented_side",
      "q_both_sides_within_library_factor", "n_features_in_set"
    )
  )
  assert(all(p$family == family), paste("Family mismatch in enrichment for", family))

  key_side <- paste(side_map$factor, side_map$oriented_side, sep = "\r")
  idx <- match(paste(p$factor, p$oriented_side, sep = "\r"), key_side)
  assert(all(!is.na(idx)), paste("Could not map all pathways to phenotypic sides for", family))
  p$factor_role <- side_map$factor_role[idx]
  p$phenotypic_side <- side_map$phenotypic_side[idx]
  p$p_assigned_class <- side_map$p_assigned_class[idx]
  p$assignment_strong <- side_map$assignment_strong[idx]
  p
}

load_observed_membership <- function(runs, primary_pathways, family) {
  libs <- unique(as.character(primary_pathways$library))
  rows <- list()
  rr <- 1L
  all_tx_features <- character()
  tx_annotation <- NULL

  for (lib in libs) {
    rds_path <- file.path(runs$enrichment, "rds", paste0("NATIVE_ENRICHMENT_", lib, ".rds"))
    assert(file.exists(rds_path), paste("Missing frozen enrichment membership RDS:", rds_path))
    obj <- readRDS(rds_path)
    assert("feature_set_matrix" %in% names(obj), paste("feature_set_matrix missing in", basename(rds_path)))
    mat <- as.matrix(obj$feature_set_matrix)
    all_tx_features <- union(all_tx_features, colnames(mat))

    pp <- primary_pathways[primary_pathways$library == lib, , drop = FALSE]
    missing_paths <- setdiff(pp$pathway, rownames(mat))
    assert(!length(missing_paths), paste("Primary pathways missing from frozen membership for", family, lib, ":", paste(missing_paths, collapse = ", ")))

    for (i in seq_len(nrow(pp))) {
      members <- colnames(mat)[mat[pp$pathway[i], , drop = TRUE] == 1]
      if (!length(members)) next
      rows[[rr]] <- data.frame(
        family = family,
        library = lib,
        pathway = pp$pathway[i],
        factor = pp$factor[i],
        factor_role = pp$factor_role[i],
        oriented_side = pp$oriented_side[i],
        phenotypic_side = pp$phenotypic_side[i],
        feature_original = as.character(members),
        formal_pathway_membership = TRUE,
        observed_in_frozen_tx_universe = TRUE,
        stringsAsFactors = FALSE
      )
      rr <- rr + 1L
    }
  }

  # MOFA_NATIVE_ENRICHMENT_OBJECT contains the complete frozen TX annotation.
  main_rds <- file.path(runs$enrichment, "rds", "MOFA_NATIVE_ENRICHMENT_OBJECT.rds")
  assert(file.exists(main_rds), paste("Missing", main_rds))
  main_obj <- readRDS(main_rds)
  assert("tx_annotation" %in% names(main_obj), paste("tx_annotation missing for", family))
  tx_annotation <- as.data.frame(main_obj$tx_annotation, stringsAsFactors = FALSE)
  assert("feature_original" %in% names(tx_annotation), paste("tx_annotation lacks feature_original for", family))
  all_model_tx_features <- unique(as.character(tx_annotation$feature_original))

  membership <- if (length(rows)) do.call(rbind, rows) else data.frame()

  list(
    membership = membership,
    tx_annotation = tx_annotation,
    all_model_tx_features = all_model_tx_features,
    observed_membership_features = unique(as.character(membership$feature_original))
  )
}

load_driver_genes <- function(runs, side_map, family) {
  x <- read_csv_checked(
    file.path(runs$crossomic, "tables", "03B_PRIMARY_PATHWAY_SIDE_DRIVER_GENES.csv"),
    required = c("library", "pathway", "factor", "oriented_side", "feature_original", "anchor_weight")
  )
  x$family <- family
  key <- paste(side_map$factor, side_map$oriented_side, sep = "\r")
  idx <- match(paste(x$factor, x$oriented_side, sep = "\r"), key)
  assert(all(!is.na(idx)), paste("Driver side map failure for", family))
  x$factor_role <- side_map$factor_role[idx]
  x$phenotypic_side <- side_map$phenotypic_side[idx]
  x
}

# ------------------------- redundancy and clustering ------------------------
pair_metrics <- function(a, b, universe_n, drivers_a = character(), drivers_b = character()) {
  a <- unique(as.character(a)); b <- unique(as.character(b))
  ia <- intersect(a, b)
  ua <- union(a, b)
  n_shared <- length(ia)
  jac <- if (length(ua)) n_shared / length(ua) else NA_real_
  ov <- if (min(length(a), length(b)) > 0L) n_shared / min(length(a), length(b)) else NA_real_
  structural <- if (is.finite(jac) && is.finite(ov)) sqrt(jac * ov) else NA_real_
  p_overlap <- if (universe_n > 0L && length(a) <= universe_n && length(b) <= universe_n) {
    stats::phyper(n_shared - 1L, length(a), universe_n - length(a), length(b), lower.tail = FALSE)
  } else {
    NA_real_
  }

  da <- unique(as.character(drivers_a)); db <- unique(as.character(drivers_b))
  di <- intersect(da, db); du <- union(da, db)
  d_jac <- if (length(du)) length(di) / length(du) else NA_real_

  list(
    n_A = length(a), n_B = length(b), n_shared_observed = n_shared,
    shared_observed_genes = if (length(ia)) paste(sort(ia), collapse = ";") else NA_character_,
    jaccard_observed = jac,
    overlap_coefficient_observed = ov,
    structural_similarity = structural,
    p_overlap_hypergeom = p_overlap,
    n_shared_drivers = length(di),
    shared_driver_genes = if (length(di)) paste(sort(di), collapse = ";") else NA_character_,
    jaccard_drivers = d_jac
  )
}

build_redundancy_table <- function(primary_pathways, membership, drivers, all_model_tx_features, family) {
  strata <- unique(primary_pathways[, c("factor", "factor_role", "phenotypic_side"), drop = FALSE])
  rows <- list(); rr <- 1L
  universe_n <- length(unique(all_model_tx_features))

  for (s in seq_len(nrow(strata))) {
    st <- strata[s, , drop = FALSE]
    pp <- primary_pathways[
      primary_pathways$factor == st$factor &
        primary_pathways$factor_role == st$factor_role &
        primary_pathways$phenotypic_side == st$phenotypic_side,
      , drop = FALSE
    ]
    if (nrow(pp) < 2L) next

    combos <- utils::combn(seq_len(nrow(pp)), 2L)
    for (j in seq_len(ncol(combos))) {
      ia <- combos[1L, j]; ib <- combos[2L, j]
      A <- pp[ia, , drop = FALSE]; B <- pp[ib, , drop = FALSE]
      ga <- membership$feature_original[membership$library == A$library & membership$pathway == A$pathway & membership$factor == A$factor & membership$phenotypic_side == A$phenotypic_side]
      gb <- membership$feature_original[membership$library == B$library & membership$pathway == B$pathway & membership$factor == B$factor & membership$phenotypic_side == B$phenotypic_side]
      da <- drivers$feature_original[drivers$library == A$library & drivers$pathway == A$pathway & drivers$factor == A$factor & drivers$phenotypic_side == A$phenotypic_side]
      db <- drivers$feature_original[drivers$library == B$library & drivers$pathway == B$pathway & drivers$factor == B$factor & drivers$phenotypic_side == B$phenotypic_side]
      pm <- pair_metrics(ga, gb, universe_n, da, db)
      rows[[rr]] <- data.frame(
        family = family,
        factor = st$factor,
        factor_role = st$factor_role,
        phenotypic_side = st$phenotypic_side,
        library_A = A$library,
        pathway_A = A$pathway,
        library_B = B$library,
        pathway_B = B$pathway,
        n_A = pm$n_A,
        n_B = pm$n_B,
        n_shared_observed = pm$n_shared_observed,
        shared_observed_genes = pm$shared_observed_genes,
        jaccard_observed = pm$jaccard_observed,
        overlap_coefficient_observed = pm$overlap_coefficient_observed,
        structural_similarity = pm$structural_similarity,
        p_overlap_hypergeom = pm$p_overlap_hypergeom,
        n_shared_drivers = pm$n_shared_drivers,
        shared_driver_genes = pm$shared_driver_genes,
        jaccard_drivers = pm$jaccard_drivers,
        same_library = identical(as.character(A$library), as.character(B$library)),
        stringsAsFactors = FALSE
      )
      rr <- rr + 1L
    }
  }

  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$q_overlap_hypergeom_within_stratum <- NA_real_
  stratum_key <- paste(out$factor, out$factor_role, out$phenotypic_side, sep = "\r")
  for (k in unique(stratum_key)) {
    idx <- which(stratum_key == k & is.finite(out$p_overlap_hypergeom))
    if (length(idx)) out$q_overlap_hypergeom_within_stratum[idx] <- stats::p.adjust(out$p_overlap_hypergeom[idx], method = "BH")
  }
  out$redundancy_edge <- is.finite(out$q_overlap_hypergeom_within_stratum) &
    out$q_overlap_hypergeom_within_stratum <= PROGRAM_OVERLAP_Q &
    is.finite(out$structural_similarity) & out$structural_similarity > 0
  out
}

program_partition_one_stratum <- function(pp, redundancy, seed) {
  # Stable node key includes library because pathway names can theoretically collide.
  node_key <- paste(pp$factor, pp$phenotypic_side, pp$library, pp$pathway, sep = "||")
  node_key <- sort(unique(node_key))

  red <- redundancy[redundancy$redundancy_edge %in% TRUE, , drop = FALSE]
  if (!nrow(red)) {
    return(data.frame(node_key = node_key, raw_cluster = seq_along(node_key), stringsAsFactors = FALSE))
  }

  edges <- data.frame(
    from = paste(red$factor, red$phenotypic_side, red$library_A, red$pathway_A, sep = "||"),
    to = paste(red$factor, red$phenotypic_side, red$library_B, red$pathway_B, sep = "||"),
    weight = red$structural_similarity,
    stringsAsFactors = FALSE
  )
  edges <- edges[order(edges$from, edges$to), , drop = FALSE]
  vertices <- data.frame(name = node_key, stringsAsFactors = FALSE)

  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
  comps <- igraph::components(g)$membership
  raw_cluster <- rep(NA_integer_, length(node_key)); names(raw_cluster) <- node_key
  next_cluster <- 1L

  for (comp in sort(unique(comps))) {
    v <- names(comps)[comps == comp]
    sg <- igraph::induced_subgraph(g, vids = v)
    if (igraph::vcount(sg) <= 2L || igraph::ecount(sg) == 0L) {
      raw_cluster[v] <- next_cluster
      next_cluster <- next_cluster + 1L
    } else {
      set.seed(seed + as.integer(comp))
      cl <- igraph::cluster_louvain(sg, weights = igraph::E(sg)$weight)
      memb <- igraph::membership(cl)
      for (m in sort(unique(memb))) {
        vv <- names(memb)[memb == m]
        raw_cluster[vv] <- next_cluster
        next_cluster <- next_cluster + 1L
      }
    }
  }

  assert(all(!is.na(raw_cluster)), "Internal clustering error: unassigned pathway node")
  data.frame(node_key = node_key, raw_cluster = as.integer(raw_cluster[node_key]), stringsAsFactors = FALSE)
}

choose_program_representatives <- function(pp, membership_primary, redundancy) {
  programs <- unique(membership_primary$program_id)
  rows <- list(); rr <- 1L

  for (pid in programs) {
    mm <- membership_primary[membership_primary$program_id == pid, , drop = FALSE]
    keys <- paste(mm$library, mm$pathway, mm$factor, mm$phenotypic_side, sep = "||")

    centrality <- rep(0, nrow(mm))
    if (nrow(mm) > 1L) {
      for (i in seq_len(nrow(mm))) {
        sims <- numeric()
        for (j in seq_len(nrow(mm))) {
          if (i == j) next
          a_lib <- mm$library[i]; a_path <- mm$pathway[i]
          b_lib <- mm$library[j]; b_path <- mm$pathway[j]
          r <- redundancy[
            redundancy$factor == mm$factor[i] &
              redundancy$phenotypic_side == mm$phenotypic_side[i] &
              (
                (redundancy$library_A == a_lib & redundancy$pathway_A == a_path & redundancy$library_B == b_lib & redundancy$pathway_B == b_path) |
                  (redundancy$library_A == b_lib & redundancy$pathway_A == b_path & redundancy$library_B == a_lib & redundancy$pathway_B == a_path)
              ),
            , drop = FALSE
          ]
          sims <- c(sims, if (nrow(r)) r$structural_similarity[1L] else 0)
        }
        centrality[i] <- mean(sims)
      }
    }

    pidx <- match(keys, paste(pp$library, pp$pathway, pp$factor, pp$phenotypic_side, sep = "||"))
    qvals <- pp$q_both_sides_within_library_factor[pidx]
    sizes <- pp$n_features_in_set[pidx]
    ord <- order(-centrality, qvals, -sizes, mm$library, mm$pathway, na.last = TRUE)
    rep_i <- ord[1L]
    best_i <- order(qvals, -sizes, mm$library, mm$pathway, na.last = TRUE)[1L]

    hallmark_idx <- which(mm$library == "HALLMARK")
    hallmark_path <- NA_character_
    if (length(hallmark_idx)) {
      hi <- hallmark_idx[order(qvals[hallmark_idx], mm$pathway[hallmark_idx], na.last = TRUE)][1L]
      hallmark_path <- mm$pathway[hi]
    }

    rows[[rr]] <- data.frame(
      program_id = pid,
      representative_library = mm$library[rep_i],
      representative_pathway = mm$pathway[rep_i],
      representative_factor = mm$factor[rep_i],
      representative_oriented_side = mm$oriented_side[rep_i],
      representative_phenotypic_side = mm$phenotypic_side[rep_i],
      representative_mean_structural_similarity = centrality[rep_i],
      best_supported_library = mm$library[best_i],
      best_supported_pathway = mm$pathway[best_i],
      best_supported_q = qvals[best_i],
      hallmark_anchor = hallmark_path,
      n_primary_pathways = nrow(mm),
      stringsAsFactors = FALSE
    )
    rr <- rr + 1L
  }
  do.call(rbind, rows)
}

build_programs <- function(primary_pathways, redundancy, family) {
  strata <- unique(primary_pathways[, c("factor", "factor_role", "phenotypic_side"), drop = FALSE])
  primary_rows <- list(); rr <- 1L

  for (s in seq_len(nrow(strata))) {
    st <- strata[s, , drop = FALSE]
    pp <- primary_pathways[
      primary_pathways$factor == st$factor &
        primary_pathways$factor_role == st$factor_role &
        primary_pathways$phenotypic_side == st$phenotypic_side,
      , drop = FALSE
    ]
    red <- redundancy[
      redundancy$factor == st$factor &
        redundancy$factor_role == st$factor_role &
        redundancy$phenotypic_side == st$phenotypic_side,
      , drop = FALSE
    ]
    part <- program_partition_one_stratum(pp, red, PROGRAM_GRAPH_SEED + s * 100L)
    pp$node_key <- paste(pp$factor, pp$phenotypic_side, pp$library, pp$pathway, sep = "||")
    pp$raw_cluster <- part$raw_cluster[match(pp$node_key, part$node_key)]
    assert(all(!is.na(pp$raw_cluster)), paste("Program partition failure for", family, st$factor, st$phenotypic_side))

    # Stable temporary program label; final label will be replaced by representative-based ID.
    pp$tmp_program <- paste(st$factor, st$factor_role, st$phenotypic_side, sprintf("C%03d", pp$raw_cluster), sep = "__")
    primary_rows[[rr]] <- pp
    rr <- rr + 1L
  }

  ptmp <- do.call(rbind, primary_rows)
  # Convert temporary clusters to stable IDs sorted by stratum + lexicographic node set.
  clusters <- split(ptmp, ptmp$tmp_program)
  cluster_meta <- lapply(clusters, function(z) {
    data.frame(
      tmp_program = z$tmp_program[1L],
      factor = z$factor[1L],
      factor_role = z$factor_role[1L],
      phenotypic_side = z$phenotypic_side[1L],
      lexical_signature = paste(sort(paste(z$library, z$pathway, sep = "||")), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  cm <- do.call(rbind, cluster_meta)
  cm <- cm[order(cm$factor_role, cm$phenotypic_side, cm$lexical_signature), , drop = FALSE]
  cm$program_index <- ave(seq_len(nrow(cm)), paste(cm$factor_role, cm$phenotypic_side, sep = "\r"), FUN = seq_along)
  cm$program_id <- paste0(
    family, "__", cm$factor_role, "__", cm$phenotypic_side, "__P", sprintf("%02d", cm$program_index)
  )
  ptmp$program_id <- cm$program_id[match(ptmp$tmp_program, cm$tmp_program)]
  ptmp$membership_type <- "PRIMARY"
  ptmp$bridge_pathway <- FALSE

  # Secondary affiliations arise only from the same statistically supported
  # redundancy edges used in the graph. No additional similarity threshold.
  sec_rows <- list(); ss <- 1L
  if (nrow(redundancy)) {
    edges <- redundancy[redundancy$redundancy_edge %in% TRUE, , drop = FALSE]
    if (nrow(edges)) {
      key_to_program <- setNames(ptmp$program_id, ptmp$node_key)
      for (i in seq_len(nrow(edges))) {
        ka <- paste(edges$factor[i], edges$phenotypic_side[i], edges$library_A[i], edges$pathway_A[i], sep = "||")
        kb <- paste(edges$factor[i], edges$phenotypic_side[i], edges$library_B[i], edges$pathway_B[i], sep = "||")
        pa <- unname(key_to_program[ka]); pb <- unname(key_to_program[kb])
        if (is.na(pa) || is.na(pb) || pa == pb) next

        za <- ptmp[ptmp$node_key == ka, , drop = FALSE]
        zb <- ptmp[ptmp$node_key == kb, , drop = FALSE]
        za$program_id <- pb; za$membership_type <- "SECONDARY_CROSS_COMMUNITY_EDGE"; za$bridge_pathway <- TRUE
        zb$program_id <- pa; zb$membership_type <- "SECONDARY_CROSS_COMMUNITY_EDGE"; zb$bridge_pathway <- TRUE
        sec_rows[[ss]] <- za; ss <- ss + 1L
        sec_rows[[ss]] <- zb; ss <- ss + 1L
      }
    }
  }

  memberships <- ptmp
  if (length(sec_rows)) {
    memberships <- rbind(memberships, do.call(rbind, sec_rows))
    memberships <- memberships[!duplicated(memberships[, c("library", "pathway", "factor", "phenotypic_side", "program_id", "membership_type")]), , drop = FALSE]
    bridge_keys <- unique(memberships$node_key[memberships$membership_type != "PRIMARY"])
    memberships$bridge_pathway <- memberships$node_key %in% bridge_keys
  }

  list(primary = ptmp, memberships = memberships)
}

# ------------------------- gene/program summaries ---------------------------
build_program_gene_membership <- function(program_memberships, observed_membership, drivers) {
  # Include primary and secondary pathway affiliations. The membership type is
  # retained so users can distinguish core versus bridge-derived membership.
  x <- merge(
    program_memberships[, c("program_id", "library", "pathway", "factor", "factor_role", "phenotypic_side", "membership_type", "bridge_pathway")],
    observed_membership[, c("library", "pathway", "factor", "factor_role", "phenotypic_side", "feature_original")],
    by = c("library", "pathway", "factor", "factor_role", "phenotypic_side"),
    all.x = TRUE,
    sort = FALSE
  )
  dkey <- paste(drivers$library, drivers$pathway, drivers$factor, drivers$phenotypic_side, drivers$feature_original, sep = "\r")
  x$is_driver_gene <- paste(x$library, x$pathway, x$factor, x$phenotypic_side, x$feature_original, sep = "\r") %in% dkey

  key <- paste(x$program_id, x$feature_original, sep = "\r")
  spl <- split(seq_len(nrow(x)), key)
  out <- lapply(spl, function(idx) {
    z <- x[idx, , drop = FALSE]
    data.frame(
      program_id = z$program_id[1L],
      factor = z$factor[1L],
      factor_role = z$factor_role[1L],
      phenotypic_side = z$phenotypic_side[1L],
      feature_original = z$feature_original[1L],
      n_pathways = length(unique(paste(z$library, z$pathway, sep = "||"))),
      pathways = collapse_unique(paste(z$library, z$pathway, sep = "||")),
      n_primary_pathways = length(unique(paste(z$library[z$membership_type == "PRIMARY"], z$pathway[z$membership_type == "PRIMARY"], sep = "||"))),
      is_driver_in_any_pathway = any(z$is_driver_gene %in% TRUE, na.rm = TRUE),
      n_driver_pathways = length(unique(paste(z$library[z$is_driver_gene %in% TRUE], z$pathway[z$is_driver_gene %in% TRUE], sep = "||"))),
      gene_bridge_across_pathways = length(unique(paste(z$library, z$pathway, sep = "||"))) > 1L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

build_tx_coverage <- function(tx_annotation, observed_membership, program_memberships, final_biomarker_matrix = NULL) {
  tx <- tx_annotation
  tx$feature_original <- as.character(tx$feature_original)
  path_counts <- table(as.character(observed_membership$feature_original))
  tx$n_primary_pathway_memberships <- as.integer(path_counts[tx$feature_original])
  tx$n_primary_pathway_memberships[is.na(tx$n_primary_pathway_memberships)] <- 0L
  tx$in_any_primary_pathway <- tx$n_primary_pathway_memberships > 0L

  pm <- program_memberships[program_memberships$membership_type == "PRIMARY", , drop = FALSE]
  pg <- merge(
    pm[, c("program_id", "library", "pathway", "factor", "phenotypic_side")],
    observed_membership[, c("library", "pathway", "factor", "phenotypic_side", "feature_original")],
    by = c("library", "pathway", "factor", "phenotypic_side"),
    all = FALSE
  )
  prog_counts <- tapply(pg$program_id, pg$feature_original, function(z) length(unique(z)))
  tx$n_primary_programs <- as.integer(prog_counts[tx$feature_original])
  tx$n_primary_programs[is.na(tx$n_primary_programs)] <- 0L

  tx$present_in_final_biomarker_matrix <- FALSE
  if (!is.null(final_biomarker_matrix) && nrow(final_biomarker_matrix) && all(c("view_code", "feature_original") %in% names(final_biomarker_matrix))) {
    ff <- unique(as.character(final_biomarker_matrix$feature_original[final_biomarker_matrix$view_code == "tx"]))
    tx$present_in_final_biomarker_matrix <- tx$feature_original %in% ff
  }
  tx
}

# ---------------------- multimodal program contextualization ----------------
build_program_multimodal_context <- function(runs, family, program_primary) {
  links <- read_csv_checked(
    file.path(runs$crossomic, "tables", "04_PATHWAY_FEATURE_CROSSOMIC_LINKS_ALL.csv"),
    required = c("family", "library", "pathway", "factor", "oriented_side", "view_code", "feature_original", "feature_model", "pathway_feature_cosine")
  )
  train <- read_csv_checked(
    file.path(runs$train, "tables", "04_PATHWAY_FEATURE_ASSOCIATIONS_OBSERVED_TRAIN.csv"),
    required = c("library", "pathway", "factor", "oriented_side", "view_code", "feature_original", "rho_spearman_train", "direction_concordant_latent_train", "q_perm_directional_within_library_view", "jackknife_expected_sign_stability")
  )
  test <- read_csv_checked(
    file.path(runs$test, "tables", "04_PATHWAY_FEATURE_ASSOCIATIONS_TEST.csv"),
    required = c("library", "pathway", "factor", "oriented_side", "view_code", "feature_original", "feature_present_in_test", "rho_test", "latent_test_direction_concordant", "train_test_same_sign", "latent_train_test_same_sign", "transport_status")
  )

  pmap <- program_primary[, c("program_id", "library", "pathway", "factor", "oriented_side", "factor_role", "phenotypic_side"), drop = FALSE]
  x <- merge(links, pmap, by = c("library", "pathway", "factor", "oriented_side"), all = FALSE, sort = FALSE)

  train_key <- paste(train$library, train$pathway, train$factor, train$oriented_side, train$view_code, train$feature_original, sep = "\r")
  test_key <- paste(test$library, test$pathway, test$factor, test$oriented_side, test$view_code, test$feature_original, sep = "\r")
  xkey <- paste(x$library, x$pathway, x$factor, x$oriented_side, x$view_code, x$feature_original, sep = "\r")
  it <- match(xkey, train_key); ie <- match(xkey, test_key)

  x$rho_spearman_train <- train$rho_spearman_train[it]
  x$direction_concordant_latent_train <- train$direction_concordant_latent_train[it]
  x$q_perm_directional_within_library_view <- train$q_perm_directional_within_library_view[it]
  x$jackknife_expected_sign_stability <- train$jackknife_expected_sign_stability[it]
  x$feature_present_in_test <- test$feature_present_in_test[ie]
  x$rho_test <- test$rho_test[ie]
  x$latent_test_direction_concordant <- test$latent_test_direction_concordant[ie]
  x$train_test_same_sign <- test$train_test_same_sign[ie]
  x$latent_train_test_same_sign <- test$latent_train_test_same_sign[ie]
  x$transport_status <- test$transport_status[ie]

  # Program-feature aggregation: no composite score and no feature selection.
  grp <- paste(x$program_id, x$view_code, x$feature_original, x$feature_model, sep = "\r")
  spl <- split(seq_len(nrow(x)), grp)
  out <- lapply(spl, function(idx) {
    z <- x[idx, , drop = FALSE]
    data.frame(
      family = family,
      program_id = z$program_id[1L],
      factor = z$factor[1L],
      factor_role = z$factor_role[1L],
      phenotypic_side = z$phenotypic_side[1L],
      view_code = z$view_code[1L],
      feature_original = z$feature_original[1L],
      feature_model = z$feature_model[1L],
      n_program_pathways_linked = length(unique(paste(z$library, z$pathway, sep = "||"))),
      linked_pathways = collapse_unique(paste(z$library, z$pathway, sep = "||")),
      median_pathway_feature_cosine = safe_median(z$pathway_feature_cosine),
      max_abs_pathway_feature_cosine = safe_max(abs(z$pathway_feature_cosine)),
      median_rho_train = safe_median(z$rho_spearman_train),
      n_train_latent_direction_concordant = n_true(z$direction_concordant_latent_train),
      n_train_q10_directional = sum(is.finite(z$q_perm_directional_within_library_view) & z$q_perm_directional_within_library_view <= 0.10, na.rm = TRUE),
      median_train_jackknife_expected_sign_stability = safe_median(z$jackknife_expected_sign_stability),
      feature_present_in_test_any = any(as_logical_safe(z$feature_present_in_test) %in% TRUE, na.rm = TRUE),
      median_rho_test = safe_median(z$rho_test),
      n_test_latent_direction_concordant = n_true(z$latent_test_direction_concordant),
      n_train_test_same_sign = n_true(z$train_test_same_sign),
      n_latent_train_test_same_sign = n_true(z$latent_train_test_same_sign),
      test_transport_statuses = collapse_unique(z$transport_status),
      evidence_interpretation = if (z$view_code[1L] == "tx") "direct_TX_membership_or_latent_context_check_membership_table" else "program_associated_not_biochemical_membership",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

build_program_view_summary <- function(program_primary, context) {
  pids <- unique(program_primary$program_id)
  views <- c("tx", "pr", "me", "cl")
  rows <- list(); rr <- 1L
  for (pid in pids) {
    for (v in views) {
      z <- context[context$program_id == pid & context$view_code == v, , drop = FALSE]
      rows[[rr]] <- data.frame(
        program_id = pid,
        view_code = v,
        n_features_with_latent_context = nrow(z),
        n_features_with_any_train_q10_support = if (nrow(z)) sum(z$n_train_q10_directional > 0L, na.rm = TRUE) else 0L,
        n_features_with_any_test_direction_preserved = if (nrow(z)) sum(grepl("DIRECTION_PRESERVED|PRESERVED", z$test_transport_statuses), na.rm = TRUE) else 0L,
        evidence_role = if (v == "tx") "TX_direct_membership_plus_latent_context" else "program_associated_context",
        stringsAsFactors = FALSE
      )
      rr <- rr + 1L
    }
  }
  do.call(rbind, rows)
}

build_pooled_representative_profiles <- function(runs, program_def) {
  pooled <- read_csv_checked(
    file.path(runs$test, "tables", "08_POOLED_PATHWAY_SCORES_DESCRIPTIVE.csv"),
    required = c("family", "dataset", "library", "pathway", "factor", "oriented_side", "sample_id", "grupo", "pathway_score")
  )
  repmap <- program_def[, c(
    "program_id", "representative_library", "representative_pathway",
    "representative_factor", "representative_oriented_side"
  ), drop = FALSE]
  names(repmap)[2:5] <- c("library", "pathway", "factor", "oriented_side")
  x <- merge(pooled, repmap, by = c("library", "pathway", "factor", "oriented_side"), all = FALSE, sort = FALSE)
  x$profile_type <- "frozen_representative_pathway_score_only"
  x$used_for_program_definition <- FALSE
  x
}

# ---------------------------- family processing -----------------------------
process_family <- function(runs, family, out_family) {
  dir.create(out_family, recursive = TRUE, showWarnings = FALSE)

  side_post <- read_csv_checked(
    file.path(runs$class_pathway, "tables", "01_FACTOR_SIDE_CLASS_POSTERIOR.csv")
  )
  side_map <- infer_factor_roles_and_sides(side_post, family)
  write_csv(side_map, file.path(out_family, "01_FACTOR_PHENOTYPIC_SIDE_MAP.csv"))

  primary <- load_primary_pathways(runs, family, side_map)
  write_csv(primary, file.path(out_family, "02_PRIMARY_PATHWAYS_WITH_PHENOTYPIC_SIDE.csv"))

  mem_obj <- load_observed_membership(runs, primary, family)
  observed <- mem_obj$membership
  drivers <- load_driver_genes(runs, side_map, family)

  dkey <- paste(drivers$library, drivers$pathway, drivers$factor, drivers$phenotypic_side, drivers$feature_original, sep = "\r")
  observed$is_side_driver <- paste(observed$library, observed$pathway, observed$factor, observed$phenotypic_side, observed$feature_original, sep = "\r") %in% dkey
  write_csv(observed, file.path(out_family, "03_PRIMARY_PATHWAY_GENE_OBSERVED_MEMBERSHIP_LONG.csv"))
  write_csv(drivers, file.path(out_family, "04_PRIMARY_PATHWAY_SIDE_DRIVER_GENES.csv"))

  final_biomarker <- read_csv_checked(file.path(runs$test, "tables", "07_BIOMARKER_EVIDENCE_MATRIX_FINAL.csv"))

  redundancy <- build_redundancy_table(primary, observed, drivers, mem_obj$all_model_tx_features, family)
  write_csv(redundancy, file.path(out_family, "05_PATHWAY_PAIRWISE_REDUNDANCY.csv"))

  progs <- build_programs(primary, redundancy, family)
  primary_members <- progs$primary
  memberships <- progs$memberships

  reps <- choose_program_representatives(primary, primary_members, redundancy)
  # Attach representatives and replace IDs only by stable already generated IDs.
  program_def <- unique(primary_members[, c("family", "factor", "factor_role", "phenotypic_side", "program_id"), drop = FALSE])
  program_def <- merge(program_def, reps, by = "program_id", all.x = TRUE, sort = FALSE)
  program_def$n_secondary_affiliations <- vapply(program_def$program_id, function(pid) {
    sum(memberships$program_id == pid & memberships$membership_type != "PRIMARY")
  }, integer(1))
  program_def$program_definition_rule <- paste0(
    "frozen_Gobserved; hypergeom_BH_q<=", format(PROGRAM_OVERLAP_Q, scientific = FALSE),
    "; edge_weight=sqrt(Jaccard*OverlapCoeff); Louvain_with_fixed_seed; no_TEST; no_PR_ME_CL"
  )

  # Add representative info to pathway membership.
  memberships <- merge(
    memberships,
    reps[, c("program_id", "representative_library", "representative_pathway")],
    by = "program_id",
    all.x = TRUE,
    sort = FALSE
  )
  memberships <- memberships[order(memberships$program_id, memberships$membership_type, memberships$library, memberships$pathway), , drop = FALSE]

  write_csv(program_def, file.path(out_family, "06_TX_PROGRAM_DEFINITION.csv"))
  write_csv(memberships, file.path(out_family, "07_TX_PROGRAM_PATHWAY_MEMBERSHIP.csv"))

  program_genes <- build_program_gene_membership(memberships, observed, drivers)
  write_csv(program_genes, file.path(out_family, "08_TX_PROGRAM_GENE_MEMBERSHIP.csv"))

  tx_cov <- build_tx_coverage(mem_obj$tx_annotation, observed, memberships, final_biomarker)
  write_csv(tx_cov, file.path(out_family, "09_TX_FEATURE_PRIMARY_PATHWAY_COVERAGE.csv"))
  write_csv(
    tx_cov[!tx_cov$in_any_primary_pathway, , drop = FALSE],
    file.path(out_family, "09B_TX_FEATURES_WITHOUT_PRIMARY_PATHWAY.csv")
  )

  context <- build_program_multimodal_context(runs, family, primary_members)
  write_csv(context, file.path(out_family, "10_PROGRAM_MULTIMODAL_FEATURE_CONTEXT.csv"))

  view_summary <- build_program_view_summary(primary_members, context)
  view_summary$family <- family
  write_csv(view_summary, file.path(out_family, "11_PROGRAM_VIEW_EVIDENCE_SUMMARY.csv"))

  pooled_rep <- build_pooled_representative_profiles(runs, program_def)
  write_csv(pooled_rep, file.path(out_family, "12_POOLED_FROZEN_PROGRAM_REPRESENTATIVE_PROFILES.csv"))

  # Run audit table.
  audit <- data.frame(
    family = family,
    bidirectional_run = runs$bidirectional,
    latent_run = runs$latent,
    enrichment_run = runs$enrichment,
    class_pathway_run = runs$class_pathway,
    crossomic_run = runs$crossomic,
    observed_train_run = runs$train,
    test_run = runs$test,
    n_primary_pathways = nrow(primary),
    n_tx_programs = nrow(program_def),
    n_tx_features_model = length(mem_obj$all_model_tx_features),
    n_tx_features_in_any_primary_pathway = length(unique(observed$feature_original)),
    n_tx_features_without_primary_pathway = sum(!tx_cov$in_any_primary_pathway),
    program_overlap_q = PROGRAM_OVERLAP_Q,
    graph_seed = PROGRAM_GRAPH_SEED,
    igraph_version = as.character(utils::packageVersion("igraph")),
    stringsAsFactors = FALSE
  )
  write_csv(audit, file.path(out_family, "00_PROGRAM_DEFINITION_AUDIT.csv"))

  list(
    family = family,
    side_map = side_map,
    primary = primary,
    observed = observed,
    drivers = drivers,
    redundancy = redundancy,
    program_def = program_def,
    memberships = memberships,
    program_genes = program_genes,
    tx_coverage = tx_cov,
    context = context,
    view_summary = view_summary,
    pooled_rep = pooled_rep,
    audit = audit
  )
}

# ---------------------------- cross-family layer ----------------------------
program_gene_sets_primary <- function(res) {
  m <- res$memberships[res$memberships$membership_type == "PRIMARY", , drop = FALSE]
  x <- merge(
    m[, c("program_id", "library", "pathway", "factor", "phenotypic_side")],
    res$observed[, c("library", "pathway", "factor", "phenotypic_side", "feature_original")],
    by = c("library", "pathway", "factor", "phenotypic_side"), all = FALSE
  )
  u <- unique(x[, c("program_id", "feature_original")])
  split(u, u$program_id)
}

build_cross_family_alignment <- function(results) {
  fams <- names(results)
  if (length(fams) < 2L) return(data.frame())
  pairs <- utils::combn(fams, 2L)
  rows <- list(); rr <- 1L

  for (j in seq_len(ncol(pairs))) {
    fa <- pairs[1L, j]; fb <- pairs[2L, j]
    A <- results[[fa]]; B <- results[[fb]]
    Ad <- A$program_def; Bd <- B$program_def
    Ag <- program_gene_sets_primary(A); Bg <- program_gene_sets_primary(B)

    for (i in seq_len(nrow(Ad))) {
      for (k in seq_len(nrow(Bd))) {
        # Cross-family alignment is phenotype-based, never mathematical sign-based.
        if (Ad$factor_role[i] != Bd$factor_role[k] || Ad$phenotypic_side[i] != Bd$phenotypic_side[k]) next
        ga <- Ag[[Ad$program_id[i]]]$feature_original
        gb <- Bg[[Bd$program_id[k]]]$feature_original
        pm <- pair_metrics(ga, gb, universe_n = length(union(ga, gb)))

        pa <- A$memberships$pathway[A$memberships$program_id == Ad$program_id[i] & A$memberships$membership_type == "PRIMARY"]
        pb <- B$memberships$pathway[B$memberships$program_id == Bd$program_id[k] & B$memberships$membership_type == "PRIMARY"]
        shared_paths <- intersect(unique(pa), unique(pb))
        path_union <- union(unique(pa), unique(pb))

        rows[[rr]] <- data.frame(
          family_A = fa,
          family_B = fb,
          factor_role = Ad$factor_role[i],
          phenotypic_side = Ad$phenotypic_side[i],
          program_A = Ad$program_id[i],
          program_B = Bd$program_id[k],
          representative_A = Ad$representative_pathway[i],
          representative_B = Bd$representative_pathway[k],
          n_genes_A = length(unique(ga)),
          n_genes_B = length(unique(gb)),
          n_shared_observed_genes = length(intersect(ga, gb)),
          jaccard_observed_genes = if (length(union(ga, gb))) length(intersect(ga, gb)) / length(union(ga, gb)) else NA_real_,
          overlap_coefficient_observed_genes = if (min(length(unique(ga)), length(unique(gb))) > 0L) length(intersect(ga, gb)) / min(length(unique(ga)), length(unique(gb))) else NA_real_,
          n_shared_pathway_names = length(shared_paths),
          shared_pathway_names = if (length(shared_paths)) paste(sort(shared_paths), collapse = ";") else NA_character_,
          jaccard_pathway_names = if (length(path_union)) length(shared_paths) / length(path_union) else NA_real_,
          stringsAsFactors = FALSE
        )
        rr <- rr + 1L
      }
    }
  }

  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)

  # Reciprocal-best is a descriptive matching aid, not ranking of analyses.
  out$best_for_A <- FALSE; out$best_for_B <- FALSE
  keyA <- paste(out$family_A, out$family_B, out$program_A, sep = "\r")
  for (ka in unique(keyA)) {
    idx <- which(keyA == ka)
    mx <- max(out$jaccard_observed_genes[idx], na.rm = TRUE)
    if (is.finite(mx)) out$best_for_A[idx[out$jaccard_observed_genes[idx] == mx]] <- TRUE
  }
  keyB <- paste(out$family_A, out$family_B, out$program_B, sep = "\r")
  for (kb in unique(keyB)) {
    idx <- which(keyB == kb)
    mx <- max(out$jaccard_observed_genes[idx], na.rm = TRUE)
    if (is.finite(mx)) out$best_for_B[idx[out$jaccard_observed_genes[idx] == mx]] <- TRUE
  }
  out$reciprocal_best_match <- out$best_for_A & out$best_for_B
  out
}

build_p2_p3_sensitivity <- function(results, alignment) {
  if (!all(c("P2", "P3") %in% names(results))) return(data.frame())
  z <- alignment[alignment$family_A == "P2" & alignment$family_B == "P3", , drop = FALSE]
  if (!nrow(z)) {
    z <- alignment[alignment$family_A == "P3" & alignment$family_B == "P2", , drop = FALSE]
  }
  if (!nrow(z)) return(data.frame())
  z$comparison_scope <- "sensitivity_to_clinical_specification_with_vs_without_FLI_not_single_variable_ablation"
  z$independent_validation <- FALSE
  z
}

# ----------------------------------- main -----------------------------------
main <- function() {
  assert(length(FAMILIES) >= 1L, "No families requested")
  assert(is.finite(PROGRAM_OVERLAP_Q) && PROGRAM_OVERLAP_Q > 0 && PROGRAM_OVERLAP_Q < 1, "BIO_PROGRAM_OVERLAP_Q must be between 0 and 1")

  selection_run <- if (nzchar(SELECTION_RUN)) canonical_path(SELECTION_RUN, TRUE) else find_latest_selection_run(UPSTREAM_DIR)
  timestamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
  run_id <- paste0(timestamp, "Z_pid", Sys.getpid())
  root <- file.path(selection_run, "BIOLOGICAL_PROGRAMS", run_id)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "cross_family"), recursive = TRUE, showWarnings = FALSE)

  message("Selection run: ", selection_run)
  message("Biological-program output: ", root)
  message("Families: ", paste(FAMILIES, collapse = ", "))

  results <- list()
  manifest_rows <- list()

  for (family in FAMILIES) {
    message("\n===== ", family, " =====")
    runs <- resolve_family_runs(selection_run, family)
    out_family <- file.path(root, family)
    res <- process_family(runs, family, out_family)
    results[[family]] <- res
    manifest_rows[[family]] <- res$audit
    message("Primary pathways: ", nrow(res$primary), " | TX programs: ", nrow(res$program_def))
  }

  manifest <- do.call(rbind, manifest_rows)
  write_csv(manifest, file.path(root, "00_BIOLOGICAL_PROGRAMS_MANIFEST.csv"))

  alignment <- build_cross_family_alignment(results)
  write_csv(alignment, file.path(root, "cross_family", "01_PROGRAM_PAIRWISE_ALIGNMENT_PHENOTYPE_BASED.csv"))

  p2p3 <- build_p2_p3_sensitivity(results, alignment)
  write_csv(p2p3, file.path(root, "cross_family", "02_P2_P3_CLINICAL_SPECIFICATION_SENSITIVITY.csv"))

  # Save an RDS contract for script 09. This is the frozen analytical definition.
  saveRDS(
    list(
      selection_run = selection_run,
      output_root = root,
      families = FAMILIES,
      parameters = list(
        program_overlap_q = PROGRAM_OVERLAP_Q,
        graph_seed = PROGRAM_GRAPH_SEED,
        primary_pathway_definition = "upstream q_both_sides_within_library_factor <= 0.10",
        clustering_membership = "frozen observed transcriptomic pathway membership",
        clustering_test = "hypergeometric overlap + BH within factor-role/phenotypic-side stratum",
        edge_weight = "sqrt(Jaccard_observed * overlap_coefficient_observed)",
        community_algorithm = "igraph Louvain with fixed seed; singleton programs allowed",
        driver_genes_used_for_clustering = FALSE,
        PR_ME_CL_used_for_TX_program_definition = FALSE,
        TEST_used_for_program_definition = FALSE,
        pooled_used_for_program_definition = FALSE,
        P2_P3_interpretation = "clinical specification sensitivity with/without FLI, not single-variable causal ablation"
      ),
      family_results = lapply(results, function(x) {
        list(
          side_map = x$side_map,
          program_def = x$program_def,
          memberships = x$memberships,
          program_genes = x$program_genes,
          tx_coverage = x$tx_coverage,
          context = x$context,
          view_summary = x$view_summary
        )
      }),
      cross_family_alignment = alignment,
      p2_p3_sensitivity = p2p3
    ),
    file.path(root, "BIOLOGICAL_PROGRAMS_OBJECT.rds"),
    version = 3
  )

  readme <- c(
    "METATEST — BIOLOGICAL PROGRAMS",
    "",
    "This run defines TX programs and attaches multimodal context from scripts 01–07.",
    "It does not create a biomarker shortlist and does not build final figures.",
    "",
    "Key reproducibility rules:",
    "- Phenotypic sides are reconstructed from the frozen class-factor posterior audit; mathematical sign is never assumed to carry a fixed phenotype.",
    "- TX programs use frozen observed formal pathway membership only.",
    "- Drivers support interpretation but do not create clusters.",
    "- PR/ME/CL are attached only after TX programs are frozen.",
    "- TEST is transportability only; pooled n=26 is descriptive only.",
    "- P2/P3 comparison is clinical-specification sensitivity with/without FLI, not a pure single-variable FLI ablation.",
    "- TX genes with no primary pathway membership are preserved in 09_TX_FEATURE_PRIMARY_PATHWAY_COVERAGE.csv.",
    "",
    paste0("PROGRAM_OVERLAP_Q=", PROGRAM_OVERLAP_Q),
    paste0("PROGRAM_GRAPH_SEED=", PROGRAM_GRAPH_SEED),
    paste0("igraph_version=", as.character(utils::packageVersion("igraph")))
  )
  writeLines(readme, file.path(root, "READ_ME_FIRST.txt"), useBytes = TRUE)

  latest_root <- file.path(selection_run, "BIOLOGICAL_PROGRAMS")
  dir.create(latest_root, recursive = TRUE, showWarnings = FALSE)
  writeLines(canonical_path(root, TRUE), file.path(latest_root, "LATEST_RUN.txt"), useBytes = TRUE)

  message("\nDONE: ", root)
  invisible(root)
}

main()

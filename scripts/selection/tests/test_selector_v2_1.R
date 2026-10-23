#!/usr/bin/env Rscript
# Prueba sintetica de la nueva logica, sin MOFA ni datos reales.
args <- commandArgs(FALSE)
script_path <- sub('^--file=', '', args[grep('^--file=', args)][1L])
script <- file.path(dirname(dirname(normalizePath(script_path))),
                    'select_best_upstream_scenario_downstream_v2_1.R')
source(script)
mk <- function(label, q, coverage, sep, sw, sm) {
  data.frame(family='P1', scenario=label, feasible=TRUE, baseline=FALSE,
    pairwise_separation_worst=sep, silhouette_worst=sw,
    silhouette_macro=sm, dispersion_balance=q, view_balance_entropy=q,
    n_features_transcriptomica=as.integer(100*coverage),
    n_features_proteomica=23L, n_features_metabolomica=as.integer(100*coverage),
    n_features_clinical=36L)
}
a <- rbind(
 mk('MEJOR',.65,.20,.70,.55,.68),
 mk('CERCANO',.64,.45,.69,.55,.67),
 mk('LEJANO_CON_MAS_FEATURES',.20,1,.19,.01,.25))
res <- rank_scenarios(a)
stopifnot(nrow(res)==3L, sum(res$quality_plateau)==2L,
  !res$quality_plateau[res$scenario=='LEJANO_CON_MAS_FEATURES'],
  res$scenario[res$selected_scenario]=='CERCANO',
  all(c('quality_gap_to_optimum','quality_relative_to_optimum')%in%names(res)))
cat('PASS: lejano excluido, cobertura opera solo dentro de la zona cercana al optimo.\n')

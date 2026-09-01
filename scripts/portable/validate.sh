validate_heatmap_signal() {
  local dir="$1" signal="$2"
  local manifest="${dir}/diagnostico_bayes_features/final_annotated_bidirectional_module_heatmap_${signal}/tables/00_final_annotated_bidirectional_module_heatmap_manifest_${signal}.csv"
  local signature="${dir}/diagnostico_bayes_features/framework_bidireccional_${signal}/tables/01_bidirectional_feature_signature.csv"
  if [[ -s "$manifest" ]]; then
    return 0
  fi
  [[ -s "$signature" ]] || return 1
  Rscript --vanilla - "$signature" <<'RS_HEAT'
args <- commandArgs(trailingOnly = TRUE)
x <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
required <- c("evidence_grade", "bidirectional_pattern")
if (length(setdiff(required, colnames(x))) > 0L) quit(status = 1L)
eligible <- x$evidence_grade %in% c(
  "A_strong", "B_moderate", "C_exploratory", "D_ambiguous_dependent"
) & !is.na(x$bidirectional_pattern) & x$bidirectional_pattern != "uncertain"
eligible[is.na(eligible)] <- FALSE
quit(status = if (any(eligible)) 1L else 0L)
RS_HEAT
}

validate_framework_dir() {
  local dir="$1" signal relative id_resolution
  for signal in reconstruido observado; do
    for relative in \
      "framework_bidireccional_${signal}/tables/01_bidirectional_feature_signature.csv" \
      "framework_bidireccional_${signal}/tables/04_bidirectional_final_table_corrected.csv" \
      "framework_bidireccional_${signal}/plots_interpretables_bidireccional/00_plot_manifest.csv"
    do
      [[ -s "${dir}/diagnostico_bayes_features/${relative}" ]] || return 1
    done
    validate_heatmap_signal "$dir" "$signal" || return 1
  done
  [[ -s "${dir}/diagnostico_bayes_features/mofa_dual_feature_differential_ranking/20_observed_all_feature_differential_ranking.csv" ]] || return 1
  [[ -s "${dir}/diagnostico_bayes_features/mofa_dual_feature_differential_ranking/21_reconstructed_all_feature_differential_ranking.csv" ]] || return 1
  [[ -s "${dir}/diagnostico_bayes_features/database_biological_annotation/tables/00_database_annotation_manifest.csv" ]] || return 1
  id_resolution="${dir}/diagnostico_bayes_features/database_biological_annotation/tables/31_tx_pr_identifier_resolution.csv"
  [[ -s "$id_resolution" ]] || return 1
  Rscript --vanilla - "$id_resolution" <<'RS_ID'
args <- commandArgs(trailingOnly = TRUE)
x <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
quit(status = if ("selected_query" %in% colnames(x)) 0L else 1L)
RS_ID
}

#!/usr/bin/env Rscript

## ============================================================
## MAPA MAESTRO PROTEÍNAS / CITOQUINAS / MIOQUINAS
## ------------------------------------------------------------
## No depende de ningún escenario.
## No lee feature_map.
## No filtra por biomarcadores finales.
##
## Salidas:
##   reference_maps/protein_master_map.tsv
##   reference_maps/protein_master_alias_map.tsv
##   reference_maps/protein_manual_map_master.csv
##   reference_maps/protein_enrichment_gene_sets.tsv
##
## Uso:
##   Rscript crear_mapa_proteico_maestro.R
##   Rscript crear_mapa_proteico_maestro.R ./reference_maps
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
OUTDIR <- ifelse(length(args) >= 1, args[1], "reference_maps")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ============================================================
## 1) Helpers
## ============================================================

norm_protein_key <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  
  x <- gsub("^pr_", "", x)
  x <- gsub("β", "b", x)
  x <- gsub("α", "a", x)
  x <- gsub("gamma", "g", x)
  x <- gsub("interferon", "ifn", x)
  x <- gsub("tumor necrosis factor", "tnf", x)
  
  x <- gsub("\\+", " plus ", x)
  x <- gsub("/", " ", x)
  x <- gsub("-", "_", x)
  x <- gsub("[()]", " ", x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  
  ## Canonicalización de variantes frecuentes
  x <- dplyr::case_when(
    x %in% c("il_1beta", "il1beta", "il_1b", "il1b") ~ "il_1b",
    x %in% c("il_1ra", "il1ra", "il_1rn", "il1rn") ~ "il_1ra",
    x %in% c("ifng", "ifn_g", "ifn_gamma", "interferon_g") ~ "ifng",
    x %in% c("tnfa", "tnf_a", "tnf_alpha") ~ "tnfa",
    x %in% c("tnfb", "tnf_b", "tnf_beta", "lymphotoxin_alpha") ~ "tnfb",
    x %in% c("vegf_a", "vegfa") ~ "vegf_a",
    x %in% c("ip_10", "ip10", "cxcl10") ~ "ip_10",
    x %in% c("mcp_1", "mcp1", "ccl2") ~ "mcp_1",
    x %in% c("fractalkine", "fraktaline", "fractaline", "cx3cl1") ~ "fractalkine",
    x %in% c("osteocrin_musclin", "osteocrin", "musclin", "ostn") ~ "osteocrin_musclin",
    x %in% c("fstl_1", "fstl1", "follistatin_like_1") ~ "fstl_1",
    TRUE ~ x
  )
  
  x
}

sanitize_feature_suffix <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- gsub("β", "b", x)
  x <- gsub("α", "a", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

make_pr_feature_model <- function(x) {
  paste0("pr_", sanitize_feature_suffix(x))
}

collapse_chr <- function(x, sep = " | ") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return("")
  paste(x, collapse = sep)
}

## ============================================================
## 2) Diccionario maestro
## ------------------------------------------------------------
## enrichment_symbol:
##   símbolo que entra al enrichment.
##
## assay_family:
##   bloque biológico general.
##
## biological_axis:
##   eje interpretativo para red feature-pathway-grupo.
## ============================================================

protein_master <- tribble(
  ~protein_id, ~assay_label, ~canonical_protein_name, ~mapped_gene_symbol, ~mapped_uniprot_id, ~assay_family, ~biological_axis, ~enrichment_include, ~enrichment_symbol, ~mapping_confidence, ~notes,
  
  "PROT_001", "Fractalkine", "C-X3-C motif chemokine ligand 1 / Fractalkine", "CX3CL1", "P78423", "chemokine", "chemotaxis_inflammation_axis", TRUE, "CX3CL1", "high_curated_uniprot_supported", "Fractalkine = CX3CL1. Incluir alias Fraktaline/Fractaline por errores de escritura.",
  
  "PROT_002", "BDNF", "brain-derived neurotrophic factor", "BDNF", "P23560", "neurotrophin_myokine_related", "neurotrophic_secreted_factor_axis", TRUE, "BDNF", "high_curated_uniprot_supported", "Factor neurotrófico secretado; puede aparecer en paneles de miokinas.",
  
  "PROT_003", "EPO", "erythropoietin", "EPO", "P01588", "growth_factor_hormone", "erythropoiesis_hypoxia_axis", TRUE, "EPO", "high_curated_uniprot_supported", "Eritropoyetina.",
  
  "PROT_004", "Osteonectin", "secreted protein acidic and cysteine rich / Osteonectin", "SPARC", "P09486", "extracellular_matrix_secreted", "extracellular_matrix_axis", TRUE, "SPARC", "high_curated_uniprot_supported", "Osteonectin = SPARC. No usar OSTEONECTIN como símbolo.",
  
  "PROT_005", "FABP3", "fatty acid binding protein 3", "FABP3", "P05413", "metabolic_myokine_related", "fatty_acid_binding_metabolism_axis", TRUE, "FABP3", "high_curated_uniprot_supported", "Proteína de unión a ácidos grasos, músculo/corazón.",
  
  "PROT_006", "FSTL-1", "follistatin-like protein 1", "FSTL1", "Q12841", "myokine_secreted_factor", "myokine_inflammation_repair_axis", TRUE, "FSTL1", "high_curated_uniprot_supported", "FSTL-1 debe mapear a FSTL1.",
  
  "PROT_007", "FGF21", "fibroblast growth factor 21", "FGF21", "Q9NSA1", "metabolic_hormone_myokine", "energy_metabolism_axis", TRUE, "FGF21", "high_curated_uniprot_supported", "Hormona metabólica asociada a eje energético.",
  
  "PROT_008", "Osteocrin_Musclin", "osteocrin / musclin", "OSTN", "P61366", "myokine_secreted_peptide", "muscle_bone_secreted_factor_axis", TRUE, "OSTN", "high_curated_uniprot_supported", "Osteocrin/Musclin = OSTN. No usar MUSCLIN como símbolo principal.",
  
  "PROT_009", "Eotaxin", "C-C motif chemokine 11 / Eotaxin-1", "CCL11", "P51671", "chemokine", "eosinophil_chemotaxis_inflammation_axis", TRUE, "CCL11", "high_curated_uniprot_supported", "Eotaxin sin número suele referirse a Eotaxin-1 = CCL11.",
  
  "PROT_010", "IFNg", "interferon gamma", "IFNG", "P01579", "cytokine", "interferon_inflammation_axis", TRUE, "IFNG", "high_curated_uniprot_supported", "IFNg / IFN-gamma = IFNG.",
  
  "PROT_011", "IL-1b", "interleukin 1 beta", "IL1B", "P01584", "interleukin", "proinflammatory_cytokine_axis", TRUE, "IL1B", "high_curated_uniprot_supported", "IL-1b / IL-1β = IL1B.",
  
  "PROT_012", "IL-1RA", "interleukin 1 receptor antagonist", "IL1RN", "P18510", "interleukin_antagonist", "antiinflammatory_cytokine_regulation_axis", TRUE, "IL1RN", "high_curated_uniprot_supported", "IL-1RA = IL1RN, no IL1A/IL1B.",
  
  "PROT_013", "IL-2", "interleukin 2", "IL2", "P60568", "interleukin", "t_cell_cytokine_axis", TRUE, "IL2", "high_curated_uniprot_supported", "IL-2 = IL2.",
  
  "PROT_014", "IL-4", "interleukin 4", "IL4", "P05112", "interleukin", "type2_immunity_axis", TRUE, "IL4", "high_curated_uniprot_supported", "IL-4 = IL4.",
  
  "PROT_015", "IL-8", "C-X-C motif chemokine ligand 8 / interleukin 8", "CXCL8", "P10145", "chemokine_interleukin", "neutrophil_chemotaxis_inflammation_axis", TRUE, "CXCL8", "high_curated_uniprot_supported", "IL-8 se debe enriquecer como CXCL8, no como IL8 si la base usa símbolo HGNC actual.",
  
  "PROT_016", "IL-10", "interleukin 10", "IL10", "P22301", "interleukin", "antiinflammatory_cytokine_axis", TRUE, "IL10", "high_curated_uniprot_supported", "IL-10 = IL10.",
  
  "PROT_017", "IL-18", "interleukin 18", "IL18", "Q14116", "interleukin", "inflammasome_related_cytokine_axis", TRUE, "IL18", "high_curated_uniprot_supported", "IL-18 = IL18.",
  
  "PROT_018", "IL-22", "interleukin 22", "IL22", "Q9GZX6", "interleukin", "epithelial_barrier_cytokine_axis", TRUE, "IL22", "high_curated_uniprot_supported", "IL-22 = IL22.",
  
  "PROT_019", "IP-10", "C-X-C motif chemokine ligand 10 / IP-10", "CXCL10", "P02778", "chemokine", "interferon_induced_chemokine_axis", TRUE, "CXCL10", "high_curated_uniprot_supported", "IP-10 = CXCL10.",
  
  "PROT_020", "MCP-1", "C-C motif chemokine 2 / monocyte chemoattractant protein 1", "CCL2", "P13500", "chemokine", "monocyte_chemotaxis_inflammation_axis", TRUE, "CCL2", "high_curated_uniprot_supported", "MCP-1 = CCL2.",
  
  "PROT_021", "TNFa", "tumor necrosis factor alpha", "TNF", "P01375", "cytokine", "tnf_inflammation_axis", TRUE, "TNF", "high_curated_uniprot_supported", "TNFa / TNF-alpha = TNF.",
  
  "PROT_022", "TNFb", "lymphotoxin alpha / tumor necrosis factor beta", "LTA", "P01374", "cytokine", "lymphotoxin_tnf_family_axis", TRUE, "LTA", "high_curated_uniprot_supported", "TNFb / TNF-beta = LTA.",
  
  "PROT_023", "VEGF-A", "vascular endothelial growth factor A", "VEGFA", "P15692", "growth_factor", "angiogenesis_axis", TRUE, "VEGFA", "high_curated_uniprot_supported", "VEGF-A = VEGFA."
) %>%
  mutate(
    database_id = mapped_uniprot_id,
    database_id_type = "UniProtKB",
    feature_model_suggested = make_pr_feature_model(assay_label),
    feature_key = norm_protein_key(assay_label),
    query_key = norm_protein_key(canonical_protein_name)
  )

## ============================================================
## 3) Alias amplios
## ============================================================

alias_raw <- tribble(
  ~protein_id, ~alias,
  
  "PROT_001", "Fractalkine",
  "PROT_001", "Fraktaline",
  "PROT_001", "Fractaline",
  "PROT_001", "CX3CL1",
  "PROT_001", "pr_Fractalkine",
  "PROT_001", "pr_Fraktaline",
  "PROT_001", "pr_CX3CL1",
  
  "PROT_002", "BDNF",
  "PROT_002", "brain-derived neurotrophic factor",
  "PROT_002", "pr_BDNF",
  
  "PROT_003", "EPO",
  "PROT_003", "erythropoietin",
  "PROT_003", "pr_EPO",
  
  "PROT_004", "Osteonectin",
  "PROT_004", "SPARC",
  "PROT_004", "secreted protein acidic and cysteine rich",
  "PROT_004", "pr_Osteonectin",
  "PROT_004", "pr_SPARC",
  
  "PROT_005", "FABP3",
  "PROT_005", "fatty acid binding protein 3",
  "PROT_005", "pr_FABP3",
  
  "PROT_006", "FSTL-1",
  "PROT_006", "FSTL1",
  "PROT_006", "FSTL_1",
  "PROT_006", "follistatin-like 1",
  "PROT_006", "follistatin like protein 1",
  "PROT_006", "pr_FSTL-1",
  "PROT_006", "pr_FSTL_1",
  "PROT_006", "pr_FSTL1",
  
  "PROT_007", "FGF21",
  "PROT_007", "fibroblast growth factor 21",
  "PROT_007", "pr_FGF21",
  
  "PROT_008", "Osteocrin_Musclin",
  "PROT_008", "Osteocrin",
  "PROT_008", "Musclin",
  "PROT_008", "OSTN",
  "PROT_008", "pr_Osteocrin_Musclin",
  "PROT_008", "pr_Osteocrin",
  "PROT_008", "pr_Musclin",
  "PROT_008", "pr_OSTN",
  
  "PROT_009", "Eotaxin",
  "PROT_009", "Eotaxin-1",
  "PROT_009", "CCL11",
  "PROT_009", "pr_Eotaxin",
  "PROT_009", "pr_Eotaxin_1",
  "PROT_009", "pr_CCL11",
  
  "PROT_010", "IFNg",
  "PROT_010", "IFN-gamma",
  "PROT_010", "IFN_gamma",
  "PROT_010", "interferon gamma",
  "PROT_010", "IFNG",
  "PROT_010", "pr_IFNg",
  "PROT_010", "pr_IFN_gamma",
  "PROT_010", "pr_IFNG",
  
  "PROT_011", "IL-1b",
  "PROT_011", "IL_1b",
  "PROT_011", "IL1b",
  "PROT_011", "IL-1beta",
  "PROT_011", "IL1B",
  "PROT_011", "interleukin 1 beta",
  "PROT_011", "pr_IL-1b",
  "PROT_011", "pr_IL_1b",
  "PROT_011", "pr_IL1B",
  
  "PROT_012", "IL-1RA",
  "PROT_012", "IL_1RA",
  "PROT_012", "IL1RA",
  "PROT_012", "IL1RN",
  "PROT_012", "interleukin 1 receptor antagonist",
  "PROT_012", "pr_IL-1RA",
  "PROT_012", "pr_IL_1RA",
  "PROT_012", "pr_IL1RN",
  
  "PROT_013", "IL-2",
  "PROT_013", "IL_2",
  "PROT_013", "IL2",
  "PROT_013", "interleukin 2",
  "PROT_013", "pr_IL-2",
  "PROT_013", "pr_IL_2",
  "PROT_013", "pr_IL2",
  
  "PROT_014", "IL-4",
  "PROT_014", "IL_4",
  "PROT_014", "IL4",
  "PROT_014", "interleukin 4",
  "PROT_014", "pr_IL-4",
  "PROT_014", "pr_IL_4",
  "PROT_014", "pr_IL4",
  
  "PROT_015", "IL-8",
  "PROT_015", "IL_8",
  "PROT_015", "IL8",
  "PROT_015", "CXCL8",
  "PROT_015", "interleukin 8",
  "PROT_015", "pr_IL-8",
  "PROT_015", "pr_IL_8",
  "PROT_015", "pr_IL8",
  "PROT_015", "pr_CXCL8",
  
  "PROT_016", "IL-10",
  "PROT_016", "IL_10",
  "PROT_016", "IL10",
  "PROT_016", "interleukin 10",
  "PROT_016", "pr_IL-10",
  "PROT_016", "pr_IL_10",
  "PROT_016", "pr_IL10",
  
  "PROT_017", "IL-18",
  "PROT_017", "IL_18",
  "PROT_017", "IL18",
  "PROT_017", "interleukin 18",
  "PROT_017", "pr_IL-18",
  "PROT_017", "pr_IL_18",
  "PROT_017", "pr_IL18",
  
  "PROT_018", "IL-22",
  "PROT_018", "IL_22",
  "PROT_018", "IL22",
  "PROT_018", "interleukin 22",
  "PROT_018", "pr_IL-22",
  "PROT_018", "pr_IL_22",
  "PROT_018", "pr_IL22",
  
  "PROT_019", "IP-10",
  "PROT_019", "IP_10",
  "PROT_019", "IP10",
  "PROT_019", "CXCL10",
  "PROT_019", "10 kDa interferon gamma-induced protein",
  "PROT_019", "pr_IP-10",
  "PROT_019", "pr_IP_10",
  "PROT_019", "pr_IP10",
  "PROT_019", "pr_CXCL10",
  
  "PROT_020", "MCP-1",
  "PROT_020", "MCP_1",
  "PROT_020", "MCP1",
  "PROT_020", "CCL2",
  "PROT_020", "monocyte chemoattractant protein 1",
  "PROT_020", "pr_MCP-1",
  "PROT_020", "pr_MCP_1",
  "PROT_020", "pr_MCP1",
  "PROT_020", "pr_CCL2",
  
  "PROT_021", "TNFa",
  "PROT_021", "TNF-a",
  "PROT_021", "TNF_alpha",
  "PROT_021", "TNF-alpha",
  "PROT_021", "TNF",
  "PROT_021", "pr_TNFa",
  "PROT_021", "pr_TNF_alpha",
  "PROT_021", "pr_TNF",
  
  "PROT_022", "TNFb",
  "PROT_022", "TNF-b",
  "PROT_022", "TNF_beta",
  "PROT_022", "TNF-beta",
  "PROT_022", "LTA",
  "PROT_022", "lymphotoxin alpha",
  "PROT_022", "pr_TNFb",
  "PROT_022", "pr_TNF_beta",
  "PROT_022", "pr_LTA",
  
  "PROT_023", "VEGF-A",
  "PROT_023", "VEGF_A",
  "PROT_023", "VEGFA",
  "PROT_023", "vascular endothelial growth factor A",
  "PROT_023", "pr_VEGF-A",
  "PROT_023", "pr_VEGF_A",
  "PROT_023", "pr_VEGFA"
)

## Añadir alias automáticos desde assay_label, canonical name y symbol.
alias_auto <- protein_master %>%
  transmute(
    protein_id,
    alias = assay_label
  ) %>%
  bind_rows(
    protein_master %>% transmute(protein_id, alias = mapped_gene_symbol),
    protein_master %>% transmute(protein_id, alias = canonical_protein_name),
    protein_master %>% transmute(protein_id, alias = feature_model_suggested),
    protein_master %>% transmute(protein_id, alias = make_pr_feature_model(mapped_gene_symbol)),
    protein_master %>% transmute(protein_id, alias = make_pr_feature_model(assay_label))
  )

protein_alias_map <- bind_rows(alias_raw, alias_auto) %>%
  mutate(
    alias = as.character(alias),
    alias_key = norm_protein_key(alias),
    feature_model_alias = if_else(
      grepl("^pr_", alias),
      alias,
      make_pr_feature_model(alias)
    )
  ) %>%
  filter(!is.na(alias), alias != "", !is.na(alias_key), alias_key != "") %>%
  left_join(protein_master, by = "protein_id") %>%
  distinct(alias_key, .keep_all = TRUE) %>%
  arrange(protein_id, alias)

## ============================================================
## 4) Manual map compatible con tu pipeline actual
## ------------------------------------------------------------
## Tu pipeline espera:
##   feature_model
##   feature_label
##   database_id
##   database_id_type
##   mapped_gene_symbol
##   mapped_uniprot_id
##   mapping_confidence
##   notes
## ============================================================

protein_manual_map_master <- protein_alias_map %>%
  transmute(
    feature_model = feature_model_alias,
    feature_label = assay_label,
    database_id = mapped_uniprot_id,
    database_id_type = "UniProtKB",
    mapped_gene_symbol = mapped_gene_symbol,
    mapped_uniprot_id = mapped_uniprot_id,
    mapping_confidence = mapping_confidence,
    notes = paste0(
      "Curated protein/cytokine/myokine map. ",
      "Assay label: ", assay_label,
      "; canonical protein: ", canonical_protein_name,
      "; enrichment symbol: ", enrichment_symbol,
      "; family: ", assay_family,
      "; axis: ", biological_axis,
      ". ", notes
    )
  ) %>%
  distinct(feature_model, .keep_all = TRUE) %>%
  arrange(feature_model)

## También guardar una versión sin prefijo, útil para joins por feature_label.
protein_manual_map_by_label <- protein_alias_map %>%
  transmute(
    feature_label = alias,
    feature_label_key = alias_key,
    assay_label,
    mapped_gene_symbol,
    mapped_uniprot_id,
    database_id = mapped_uniprot_id,
    database_id_type = "UniProtKB",
    assay_family,
    biological_axis,
    enrichment_symbol,
    mapping_confidence,
    notes
  ) %>%
  distinct(feature_label_key, .keep_all = TRUE) %>%
  arrange(assay_label, feature_label)

## ============================================================
## 5) Query sets biológicos opcionales para enrichment por bloques
## ============================================================

protein_enrichment_gene_sets <- protein_master %>%
  filter(enrichment_include) %>%
  select(
    protein_id,
    assay_label,
    mapped_gene_symbol,
    enrichment_symbol,
    assay_family,
    biological_axis
  ) %>%
  mutate(
    gene_symbol = enrichment_symbol
  ) %>%
  select(
    assay_family,
    biological_axis,
    protein_id,
    assay_label,
    gene_symbol,
    mapped_gene_symbol
  ) %>%
  arrange(assay_family, biological_axis, gene_symbol)

## ============================================================
## 6) Auditoría básica
## ============================================================

dup_symbols <- protein_master %>%
  count(mapped_gene_symbol) %>%
  filter(n > 1)

if (nrow(dup_symbols) > 0) {
  warning(
    "Hay símbolos duplicados en protein_master: ",
    paste(dup_symbols$mapped_gene_symbol, collapse = ", ")
  )
}

missing_uniprot <- protein_master %>%
  filter(is.na(mapped_uniprot_id) | mapped_uniprot_id == "")

if (nrow(missing_uniprot) > 0) {
  warning(
    "Hay proteínas sin UniProt: ",
    paste(missing_uniprot$assay_label, collapse = ", ")
  )
}

## ============================================================
## 7) Guardar
## ============================================================

write_tsv(
  protein_master,
  file.path(OUTDIR, "protein_master_map.tsv"),
  na = ""
)

write_tsv(
  protein_alias_map,
  file.path(OUTDIR, "protein_master_alias_map.tsv"),
  na = ""
)

write_csv(
  protein_manual_map_master,
  file.path(OUTDIR, "protein_manual_map_master.csv"),
  na = ""
)

write_csv(
  protein_manual_map_by_label,
  file.path(OUTDIR, "protein_manual_map_by_label.csv"),
  na = ""
)

write_tsv(
  protein_enrichment_gene_sets,
  file.path(OUTDIR, "protein_enrichment_gene_sets.tsv"),
  na = ""
)

cat("\n============================================================\n")
cat("MAPA MAESTRO PROTEICO GLOBAL CREADO\n")
cat("No depende de escenario.\n")
cat("Outdir:", normalizePath(OUTDIR, mustWork = FALSE), "\n")
cat("Archivos:\n")
cat(" - protein_master_map.tsv\n")
cat(" - protein_master_alias_map.tsv\n")
cat(" - protein_manual_map_master.csv\n")
cat(" - protein_manual_map_by_label.csv\n")
cat(" - protein_enrichment_gene_sets.tsv\n")
cat("============================================================\n\n")

print(
  protein_master %>%
    select(
      protein_id,
      assay_label,
      canonical_protein_name,
      mapped_gene_symbol,
      mapped_uniprot_id,
      assay_family,
      biological_axis
    ),
  n = Inf
)
#!/usr/bin/env Rscript

## ============================================================
## MAPA CLÍNICO MAESTRO GLOBAL
## ------------------------------------------------------------
## No depende de ningún escenario.
## No lee feature_map.
## No filtra por biomarcadores finales.
##
## Salidas:
##   reference_maps/clinical_master_map.tsv
##   reference_maps/clinical_master_alias_map.tsv
##   reference_maps/clinical_manual_map_master.csv
##
## Uso:
##   Rscript crear_mapa_clinico_maestro.R
##   Rscript crear_mapa_clinico_maestro.R ./reference_maps
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
OUTDIR <- ifelse(length(args) >= 1, args[1], "reference_maps")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ============================================================
## 1) Normalizador de claves clínicas
## ============================================================

norm_clinical_key <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  
  x <- gsub("^cl_", "", x)
  x <- gsub("%", " percent ", x)
  x <- gsub("\\+", " plus ", x)
  x <- gsub("/", " ratio ", x)
  x <- gsub("-", "_", x)
  x <- gsub("[()]", " ", x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  
  x
}

make_cl_feature_model <- function(feature_label) {
  paste0("cl_", gsub("[^A-Za-z0-9]+", "_", feature_label)) |>
    gsub("_+", "_", x = _) |>
    gsub("_$", "", x = _)
}

## ============================================================
## 2) Diccionario maestro
## ------------------------------------------------------------
## ontology_policy:
##
## exact_ols_allowed:
##   variable clínica simple; puede buscarse en OLS usando query inglesa.
##
## curated_index_no_ols:
##   índice derivado; NO forzar OLS.
##
## curated_ratio_no_ols:
##   ratio derivado; NO forzar OLS.
##
## curated_body_composition_no_ols:
##   variable DEXA/composición corporal; puede tener concepto relacionado,
##   pero la variable exacta suele ser derivada.
## ============================================================

clinical_master <- tribble(
  ~clinical_id, ~feature_label_es, ~canonical_label_en, ~clinical_query_en, ~semantic_category, ~biological_axis, ~preferred_databases, ~ontology_policy, ~expected_keywords, ~notes,
  
  "CLIN_001", "Edad", "age", "age", "demographics", "age_demographic_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "age", "Edad cronológica.",
  
  "CLIN_002", "GGT", "gamma-glutamyl transferase measurement", "gamma-glutamyl transferase measurement", "liver_enzyme", "hepatic_metabolism_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "gamma|glutamyl|transferase|ggt", "Marcador bioquímico hepático.",
  
  "CLIN_003", "Glu", "blood glucose measurement", "blood glucose measurement", "glucose_metabolism", "glucose_insulin_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "glucose|blood glucose", "Glucosa sanguínea.",
  
  "CLIN_004", "Col total", "total cholesterol measurement", "total cholesterol measurement", "lipid_metabolism", "lipid_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "cholesterol|total", "Colesterol total.",
  
  "CLIN_005", "LDL", "low-density lipoprotein cholesterol measurement", "low-density lipoprotein cholesterol measurement", "lipid_metabolism", "lipid_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "ldl|low density|cholesterol|lipoprotein", "Colesterol LDL.",
  
  "CLIN_006", "HDL", "high-density lipoprotein cholesterol measurement", "high-density lipoprotein cholesterol measurement", "lipid_metabolism", "lipid_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "hdl|high density|cholesterol|lipoprotein", "Colesterol HDL.",
  
  "CLIN_007", "TG", "triglyceride measurement", "triglyceride measurement", "lipid_metabolism", "lipid_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "triglyceride|triglycerides", "Triglicéridos.",
  
  "CLIN_008", "PCR", "C-reactive protein measurement", "C-reactive protein measurement", "inflammation", "inflammatory_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "c reactive|crp|protein", "PCR en español = proteína C reactiva, no polymerase chain reaction.",
  
  "CLIN_009", "Insulina", "insulin measurement", "insulin measurement", "glucose_insulin_metabolism", "glucose_insulin_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "insulin", "Insulina circulante.",
  
  "CLIN_010", "HOMA-IR", "homeostatic model assessment of insulin resistance", "homeostatic model assessment of insulin resistance", "derived_index", "glucose_insulin_axis", "Manual curated concept", "curated_index_no_ols", "homeostatic|insulin|resistance|homa", "Índice derivado de resistencia a la insulina.",
  
  "CLIN_011", "QUICKI", "quantitative insulin sensitivity check index", "quantitative insulin sensitivity check index", "derived_index", "glucose_insulin_axis", "Manual curated concept", "curated_index_no_ols", "quicki|quantitative|insulin|sensitivity", "Índice derivado de sensibilidad a insulina.",
  
  "CLIN_012", "TG-INDEX", "triglyceride glucose index", "triglyceride glucose index", "derived_index", "glucose_lipid_axis", "Manual curated concept", "curated_index_no_ols", "triglyceride|glucose|index|tyg", "Probablemente TyG index; confirmar que TG-INDEX significa triglyceride-glucose index.",
  
  "CLIN_013", "FLI", "fatty liver index", "fatty liver index", "derived_index", "hepatic_steatosis_axis", "Manual curated concept", "curated_index_no_ols", "fatty|liver|index|steatosis", "FLI = fatty liver index; nunca buscar FLI directo.",
  
  "CLIN_014", "TG/HDL", "triglyceride to HDL cholesterol ratio", "triglyceride to HDL cholesterol ratio", "derived_ratio", "lipid_insulin_resistance_axis", "Manual curated concept", "curated_ratio_no_ols", "triglyceride|hdl|ratio|cholesterol", "Ratio lipídico derivado.",
  
  "CLIN_015", "Peso", "body weight", "body weight", "anthropometry", "adiposity_axis", "EFO|HPO|NCIT|MeSH", "exact_ols_allowed", "body weight|weight", "Peso corporal; evitar falsos positivos como low-molecular-weight proteinuria.",
  
  "CLIN_016", "Altura", "body height", "body height", "anthropometry", "body_size_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "height|body height|stature", "Altura corporal.",
  
  "CLIN_017", "IMC", "body mass index", "body mass index", "anthropometry", "adiposity_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "body mass index|bmi", "IMC = body mass index; nunca buscar IMC directo.",
  
  "CLIN_018", "Ratio cintura_0M", "waist circumference related ratio", "waist circumference", "anthropometry", "central_adiposity_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "waist|circumference", "Confirmar si es perímetro cintura o ratio derivado.",
  
  "CLIN_019", "Ratio-cc", "waist-to-hip ratio", "waist-to-hip ratio", "anthropometry", "central_adiposity_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "waist|hip|ratio", "Ratio cintura-cadera.",
  
  "CLIN_020", "%BF-DEXA", "body fat percentage by dual-energy X-ray absorptiometry", "body fat percentage", "body_composition", "adiposity_axis", "Manual curated concept|EFO|NCIT", "curated_body_composition_no_ols", "body fat|fat percentage|adipose", "Porcentaje de grasa corporal medido por DEXA.",
  
  "CLIN_021", "CUN-BAE", "Clínica Universidad de Navarra body adiposity estimator", "body adiposity estimator", "derived_index", "adiposity_axis", "Manual curated concept", "curated_index_no_ols", "body|adiposity|estimator|cun bae", "CUN-BAE es un estimador de adiposidad, no un taxón.",
  
  "CLIN_022", "VAT", "visceral adipose tissue mass", "visceral adipose tissue", "body_composition", "visceral_adiposity_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "visceral|adipose|tissue", "VAT = visceral adipose tissue.",
  
  "CLIN_023", "VAT (%)", "visceral adipose tissue percentage", "visceral adipose tissue", "body_composition", "visceral_adiposity_axis", "Manual curated concept|EFO|NCIT", "curated_body_composition_no_ols", "visceral|adipose|tissue|percentage", "Porcentaje de VAT.",
  
  "CLIN_024", "Masa_grasa_0M", "fat mass", "fat mass", "body_composition", "adiposity_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "fat mass|adipose", "Masa grasa basal.",
  
  "CLIN_025", "%VAT_grasa_Total_0M", "visceral fat to total fat percentage", "visceral fat to total fat ratio", "derived_ratio", "visceral_adiposity_axis", "Manual curated concept", "curated_ratio_no_ols", "visceral|fat|total|ratio", "Ratio derivado VAT/grasa total.",
  
  "CLIN_026", "kg Ginoide", "gynoid fat mass", "gynoid fat mass", "body_composition", "fat_distribution_axis", "Manual curated concept|EFO|NCIT", "curated_body_composition_no_ols", "gynoid|fat|mass", "Masa grasa ginoide.",
  
  "CLIN_027", "% Ginoide", "gynoid fat percentage", "gynoid fat percentage", "body_composition", "fat_distribution_axis", "Manual curated concept", "curated_body_composition_no_ols", "gynoid|fat|percentage", "Porcentaje ginoide.",
  
  "CLIN_028", "% Ginoide_grasa_Total", "gynoid fat to total fat percentage", "gynoid fat to total fat ratio", "derived_ratio", "fat_distribution_axis", "Manual curated concept", "curated_ratio_no_ols", "gynoid|total|fat|ratio", "Ratio grasa ginoide/grasa total.",
  
  "CLIN_029", "kg Androide", "android fat mass", "android fat mass", "body_composition", "fat_distribution_axis", "Manual curated concept|EFO|NCIT", "curated_body_composition_no_ols", "android|fat|mass", "Masa grasa androide.",
  
  "CLIN_030", "% Androide", "android fat percentage", "android fat percentage", "body_composition", "fat_distribution_axis", "Manual curated concept", "curated_body_composition_no_ols", "android|fat|percentage", "Porcentaje androide.",
  
  "CLIN_031", "% Androide_grasa_Total", "android fat to total fat percentage", "android fat to total fat ratio", "derived_ratio", "fat_distribution_axis", "Manual curated concept", "curated_ratio_no_ols", "android|total|fat|ratio", "Ratio grasa androide/grasa total.",
  
  "CLIN_032", "Magra (%)", "lean mass percentage", "lean body mass", "body_composition", "lean_mass_axis", "Manual curated concept|EFO|NCIT", "curated_body_composition_no_ols", "lean|mass|body|percentage", "Porcentaje de masa magra.",
  
  "CLIN_033", "Magra (kg)", "lean mass", "lean body mass", "body_composition", "lean_mass_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "lean|body|mass", "Masa magra.",
  
  "CLIN_034", "Ratio and/gin", "android to gynoid fat ratio", "android to gynoid fat ratio", "derived_ratio", "fat_distribution_axis", "Manual curated concept", "curated_ratio_no_ols", "android|gynoid|ratio", "Ratio androide/ginoide.",
  
  "CLIN_035", "FFM- 0M (masa magra total + masa ósea total)", "fat-free mass", "fat-free mass", "body_composition", "lean_mass_axis", "EFO|NCIT|MeSH", "exact_ols_allowed", "fat free|lean|mass", "FFM = fat-free mass.",
  
  "CLIN_036", "Sis", "systolic blood pressure", "systolic blood pressure", "blood_pressure", "blood_pressure_axis", "EFO|HPO|NCIT|MeSH", "exact_ols_allowed", "systolic|blood pressure", "Presión arterial sistólica.",
  
  "CLIN_037", "Dias", "diastolic blood pressure", "diastolic blood pressure", "blood_pressure", "blood_pressure_axis", "EFO|HPO|NCIT|MeSH", "exact_ols_allowed", "diastolic|blood pressure", "Presión arterial diastólica."
) %>%
  mutate(
    feature_model_suggested = make_cl_feature_model(feature_label_es),
    feature_key = norm_clinical_key(feature_label_es),
    query_key = norm_clinical_key(clinical_query_en),
    database_query_type = "clinical_variable_curated",
    ontology_source = if_else(
      ontology_policy %in% c(
        "curated_index_no_ols",
        "curated_ratio_no_ols",
        "curated_body_composition_no_ols"
      ),
      "curated_clinical_concept",
      ""
    ),
    ontology_id = "",
    ontology_label = canonical_label_en,
    ontology_iri = "",
    mapping_confidence = case_when(
      ontology_policy == "exact_ols_allowed" ~ "high_curated_query_ready_for_ols_validation",
      ontology_policy == "curated_index_no_ols" ~ "high_curated_clinical_index_no_exact_ontology",
      ontology_policy == "curated_ratio_no_ols" ~ "high_curated_clinical_ratio_no_exact_ontology",
      ontology_policy == "curated_body_composition_no_ols" ~ "high_curated_body_composition_no_exact_ontology",
      TRUE ~ "requires_manual_review"
    )
  )

## ============================================================
## 3) Alias globales
## ============================================================

alias_tbl <- tribble(
  ~clinical_id, ~alias,
  
  "CLIN_001", "Edad",
  "CLIN_001", "cl_Edad",
  "CLIN_001", "age",
  
  "CLIN_002", "GGT",
  "CLIN_002", "cl_GGT",
  "CLIN_002", "gamma glutamyl transferase",
  "CLIN_002", "gamma-glutamyl transferase",
  
  "CLIN_003", "Glu",
  "CLIN_003", "cl_Glu",
  "CLIN_003", "glucosa",
  "CLIN_003", "glucose",
  
  "CLIN_004", "Col total",
  "CLIN_004", "Col_total",
  "CLIN_004", "cl_Col_total",
  "CLIN_004", "colesterol total",
  "CLIN_004", "total cholesterol",
  
  "CLIN_005", "LDL",
  "CLIN_005", "cl_LDL",
  "CLIN_005", "ldl cholesterol",
  
  "CLIN_006", "HDL",
  "CLIN_006", "cl_HDL",
  "CLIN_006", "hdl cholesterol",
  
  "CLIN_007", "TG",
  "CLIN_007", "cl_TG",
  "CLIN_007", "trigliceridos",
  "CLIN_007", "triglicéridos",
  "CLIN_007", "triglycerides",
  
  "CLIN_008", "PCR",
  "CLIN_008", "cl_PCR",
  "CLIN_008", "CRP",
  "CLIN_008", "proteina c reactiva",
  "CLIN_008", "proteína c reactiva",
  "CLIN_008", "c reactive protein",
  
  "CLIN_009", "Insulina",
  "CLIN_009", "cl_Insulina",
  "CLIN_009", "insulin",
  
  "CLIN_010", "HOMA-IR",
  "CLIN_010", "HOMA_IR",
  "CLIN_010", "cl_HOMA_IR",
  "CLIN_010", "homeostatic model assessment insulin resistance",
  
  "CLIN_011", "QUICKI",
  "CLIN_011", "cl_QUICKI",
  "CLIN_011", "quantitative insulin sensitivity check index",
  
  "CLIN_012", "TG-INDEX",
  "CLIN_012", "TG_INDEX",
  "CLIN_012", "cl_TG_INDEX",
  "CLIN_012", "TyG",
  "CLIN_012", "triglyceride glucose index",
  
  "CLIN_013", "FLI",
  "CLIN_013", "cl_FLI",
  "CLIN_013", "fatty liver index",
  
  "CLIN_014", "TG/HDL",
  "CLIN_014", "TG_HDL",
  "CLIN_014", "cl_TG_HDL",
  "CLIN_014", "triglyceride hdl ratio",
  
  "CLIN_015", "Peso",
  "CLIN_015", "cl_Peso",
  "CLIN_015", "weight",
  "CLIN_015", "body weight",
  
  "CLIN_016", "Altura",
  "CLIN_016", "cl_Altura",
  "CLIN_016", "height",
  "CLIN_016", "body height",
  
  "CLIN_017", "IMC",
  "CLIN_017", "cl_IMC",
  "CLIN_017", "BMI",
  "CLIN_017", "body mass index",
  
  "CLIN_018", "Ratio cintura_0M",
  "CLIN_018", "Ratio_cintura_0M",
  "CLIN_018", "cl_Ratio_cintura_0M",
  "CLIN_018", "waist circumference",
  "CLIN_018", "waist ratio",
  
  "CLIN_019", "Ratio-cc",
  "CLIN_019", "Ratio_cc",
  "CLIN_019", "cl_Ratio_cc",
  "CLIN_019", "waist hip ratio",
  "CLIN_019", "waist-to-hip ratio",
  
  "CLIN_020", "%BF-DEXA",
  "CLIN_020", "BF_DEXA",
  "CLIN_020", "cl_BF_DEXA",
  "CLIN_020", "body fat percentage dexa",
  
  "CLIN_021", "CUN-BAE",
  "CLIN_021", "CUN_BAE",
  "CLIN_021", "cl_CUN_BAE",
  "CLIN_021", "Clinica Universidad de Navarra body adiposity estimator",
  "CLIN_021", "Clínica Universidad de Navarra body adiposity estimator",
  
  "CLIN_022", "VAT",
  "CLIN_022", "cl_VAT",
  "CLIN_022", "visceral adipose tissue",
  
  "CLIN_023", "VAT (%)",
  "CLIN_023", "VAT_percent",
  "CLIN_023", "cl_VAT_percent",
  "CLIN_023", "visceral adipose tissue percentage",
  
  "CLIN_024", "Masa_grasa_0M",
  "CLIN_024", "cl_Masa_grasa_0M",
  "CLIN_024", "Masa grasa",
  "CLIN_024", "fat mass",
  
  "CLIN_025", "%VAT_grasa_Total_0M",
  "CLIN_025", "VAT_grasa_Total_0M",
  "CLIN_025", "cl_VAT_grasa_Total_0M",
  "CLIN_025", "visceral fat to total fat percentage",
  
  "CLIN_026", "kg Ginoide",
  "CLIN_026", "kg_Ginoide",
  "CLIN_026", "cl_kg_Ginoide",
  "CLIN_026", "gynoid fat mass",
  
  "CLIN_027", "% Ginoide",
  "CLIN_027", "percent_Ginoide",
  "CLIN_027", "cl_percent_Ginoide",
  "CLIN_027", "gynoid fat percentage",
  
  "CLIN_028", "% Ginoide_grasa_Total",
  "CLIN_028", "percent_Ginoide_grasa_Total",
  "CLIN_028", "cl_percent_Ginoide_grasa_Total",
  "CLIN_028", "gynoid fat to total fat percentage",
  
  "CLIN_029", "kg Androide",
  "CLIN_029", "kg_Androide",
  "CLIN_029", "cl_kg_Androide",
  "CLIN_029", "android fat mass",
  
  "CLIN_030", "% Androide",
  "CLIN_030", "percent_Androide",
  "CLIN_030", "cl_percent_Androide",
  "CLIN_030", "android fat percentage",
  
  "CLIN_031", "% Androide_grasa_Total",
  "CLIN_031", "percent_Androide_grasa_Total",
  "CLIN_031", "cl_percent_Androide_grasa_Total",
  "CLIN_031", "android fat to total fat percentage",
  
  "CLIN_032", "Magra (%)",
  "CLIN_032", "Magra_percent",
  "CLIN_032", "cl_Magra_percent",
  "CLIN_032", "lean mass percentage",
  
  "CLIN_033", "Magra (kg)",
  "CLIN_033", "Magra_kg",
  "CLIN_033", "cl_Magra_kg",
  "CLIN_033", "lean mass",
  
  "CLIN_034", "Ratio and/gin",
  "CLIN_034", "Ratio_and_gin",
  "CLIN_034", "cl_Ratio_and_gin",
  "CLIN_034", "android gynoid ratio",
  "CLIN_034", "android to gynoid fat ratio",
  
  "CLIN_035", "FFM- 0M (masa magra total + masa ósea total)",
  "CLIN_035", "FFM_0M",
  "CLIN_035", "cl_FFM_0M",
  "CLIN_035", "fat free mass",
  "CLIN_035", "fat-free mass",
  
  "CLIN_036", "Sis",
  "CLIN_036", "cl_Sis",
  "CLIN_036", "systolic blood pressure",
  
  "CLIN_037", "Dias",
  "CLIN_037", "cl_Dias",
  "CLIN_037", "diastolic blood pressure"
) %>%
  mutate(
    alias_key = norm_clinical_key(alias)
  ) %>%
  left_join(
    clinical_master,
    by = "clinical_id"
  ) %>%
  distinct(alias_key, .keep_all = TRUE) %>%
  arrange(clinical_id, alias)

## ============================================================
## 4) Archivo tipo manual_map compatible con anotación actual
## ============================================================

clinical_manual_map_master <- alias_tbl %>%
  filter(grepl("^cl_", alias) | alias == feature_label_es) %>%
  transmute(
    feature_model = if_else(grepl("^cl_", alias), alias, feature_model_suggested),
    feature_label = feature_label_es,
    clinical_query = clinical_query_en,
    ontology_source,
    ontology_id,
    ontology_label,
    semantic_category,
    biological_axis,
    mapping_confidence,
    notes
  ) %>%
  distinct(feature_model, .keep_all = TRUE) %>%
  arrange(feature_model)

## ============================================================
## 5) Guardar
## ============================================================

readr::write_tsv(
  clinical_master,
  file.path(OUTDIR, "clinical_master_map.tsv"),
  na = ""
)

readr::write_tsv(
  alias_tbl,
  file.path(OUTDIR, "clinical_master_alias_map.tsv"),
  na = ""
)

readr::write_csv(
  clinical_manual_map_master,
  file.path(OUTDIR, "clinical_manual_map_master.csv"),
  na = ""
)

cat("\n============================================================\n")
cat("MAPA CLÍNICO MAESTRO GLOBAL CREADO\n")
cat("No depende de escenario.\n")
cat("Outdir:", normalizePath(OUTDIR, mustWork = FALSE), "\n")
cat("Archivos:\n")
cat(" - clinical_master_map.tsv\n")
cat(" - clinical_master_alias_map.tsv\n")
cat(" - clinical_manual_map_master.csv\n")
cat("============================================================\n\n")

print(
  clinical_master %>%
    select(
      clinical_id,
      feature_label_es,
      canonical_label_en,
      clinical_query_en,
      semantic_category,
      biological_axis,
      ontology_policy
    ),
  n = Inf
)
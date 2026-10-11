# Orden oficial del descubrimiento de biomarcadores

Ejecutar con `bash scripts/orchestration/run_metatest_vigente.sh --phase biomarkers`.
Todos los scripts R incluidos aquí vienen de `scripts.zip`; el SGL se conservó
byte por byte pero recibió el prefijo **16**.

| Etapa | Archivo | Requiere |
|:--|:--|:--|
| 01 | `01_run_bidirectional_framework_selected_scenarios.R` | selección de escenarios + MOFA/RFM |
| 02 | `02_run_latent_supervised_audit.R` | 01 |
| 03 | `03_run_mofa_native_enrichment_part1.R` | 02 + msigdbr |
| 04 | `04_run_class_pathway_bridge.R` | 02 y 03 |
| 05 | `05_run_crossomic_bridge.R` | 04 |
| 06 | `06_run_observed_train_corroboration.R` | 05 |
| 07 | `07_run_test_transportability_audit.R` | 06 |
| 08 | `08_define_biological_programs.R` | 03–07 |
| 09 | `09_build_biological_synthesis_figures.R` | 08 |
| 10 | `10_build_final_predictive_panel.R` | 09 |
| 16 | `16_sparse_group_lasso_hierarchical.R` | 10 (antes `SGL_FINAL_FIXED_V3.R`) |
| 17 | `17_evidenceA_plsda_pca_validation_SCATTER_INTERSECTIONS_LOADINGS_FINAL.R` | 10 + Evidence A |
| 18 | `18_build_final_consensus_panel.R` | 17; 16 como anotación SGL |

Etapa 19_console: inspección terminal opcional. Etapa 19_SCREEN_PLOTS: inspección interactiva.
Etapas 20, 21 y 22: scripts RStudio de lectura/presentación, ejecutar manualmente
cuando el entorno gráfico esté disponible. No hacen parte de la ejecución headless.

En cada nuevo estudio, los resultados y logs se guardan fuera de Git. El
runner fija la variable `SELECTION_RUN`, pero los propios scripts deben poder
resolver las rutas internas (algunas salidas antiguas contienen rutas absolutas).

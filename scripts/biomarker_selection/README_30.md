# METATEST2 — Script 30 v1.4.1

## Diseño (9 modelos)

| Familia | MOFA_OTSU | CORE18_47 | FINAL28 |
|---|---:|---:|---:|
| P1 TODAS | 157 candidatos antes de controles TRAIN | 47 | 88 |
| P2 CON_FLI | 101 candidatos antes de controles TRAIN | 47 | 59 |
| P3 SIN_FLI | 145 candidatos antes de controles TRAIN | 47 | 67 |

- **MOFA_OTSU**: usa todas las características del modelo MOFA ganador con norma L2 de los pesos por encima de un umbral Otsu **estimado por familia y vista sobre W²**. No se exige Evidence-A ni votos predictivos. El script puede excluir variables sin varianza TRAIN o ausentes en la matriz de entrada, dejando auditoría de elegibilidad. La selección por pesos MOFA **no es supervisada por la clase**, y el umbral no es un test de significación.
- **CORE18_47**: el consenso de 47 del script 18; no se filtra de nuevo por los votos del script 28.
- **FINAL28**: todas las características finales de cada familia del script 28 (P1=88, P2=59, P3=67), sin separarlas en CORE/SHARED/SPECIFIC. El CORE de 33 del 28 queda solo como referencia descriptiva, no genera ajuste.

## Reutilización sin volver a entrenar

La versión 1.3.0 del script 30 ya ajustó 15 modelos correctamente. La versión 1.4 busca una corrida previa en el mismo directorio `SCRIPT30_LOGISTICA_MOFA` que tenga `PASS` para los 6 modelos `P1/P2/P3 × ALL/CORE18_47`, sus RDS y tablas. Reutiliza:

- `ALL` -> nuevo nombre `FINAL28` (3 ajustes).
- `CORE18_47` -> mismo nombre (3 ajustes).

Valida identidad de ruta de la exportación 28, origen del consenso 18, predictor keys y su orden, ruta del `ready_for_modeling.rds`, cantidad de muestras, número de predictores y `lambda` antes de copiar los modelos RDS y reutilizar las predicciones, coeficientes y curvas. Si falta alguna referencia, **se detiene** en vez de recalcular a escondidas.

Solo se ajustan tres modelos **MOFA_OTSU** nuevos. `MT2_30_REFIT_ALL=1` autoriza recalcular todos los 9 si fuese estrictamente necesario; evitarlo en esta corrida.

## Comandos desde raíz del repositorio

```bash
cd "$HOME/Documents/PhD/metatest/metatest"

# Smoke de UN modelo nuevo (P1, MOFA_OTSU); no reutiliza los otros.
set -o pipefail
MT2_30_SMOKE=1 Rscript --vanilla \
  scripts/biomarker_selection/30_logistica_multinomial_paneles_y_pesos_MOFA.R \
  2>&1 | tee SCRIPT30_SMOKE_v140.log
printf 'EXIT_CODE=%s\n' "${PIPESTATUS[0]}"

# Ejecución completa: 3 nuevos + 6 reutilizados.
Rscript --vanilla \
  scripts/biomarker_selection/30_logistica_multinomial_paneles_y_pesos_MOFA.R \
  2>&1 | tee SCRIPT30_v140.log
printf 'EXIT_CODE=%s\n' "${PIPESTATUS[0]}"
```

Los outputs aparecen en una nueva subcarpeta con fecha y PID, en
`results_cedia/.../seleccion_mejor_escenario_train/.../SCRIPT30_LOGISTICA_MOFA/`.

## Salidas clave

- `tables/10_MODEL_STATUS.csv`: 9 filas esperadas; `PASS` (ajuste nuevo), `PASS_REUSED` (reutilizado).
- `tables/12_MODEL_HYPERPARAMETERS.csv`, `13_MODEL_METRICS.csv`, `14_PER_CLASS_METRICS.csv`, `15_PREDICTIONS_PER_SAMPLE.csv`, `19_CONFUSION_LONG.csv`.
- `tables/23_ROC_CURVE_POINTS_OVR.csv`, `24_PR_CURVE_POINTS_OVR.csv`, `25_ROC_PR_OVR_BY_CLASS.csv`, `26_ROC_PR_MACRO_BY_MODEL.csv`.
- `tables/08_MOFA_OTSU_PANEL_MEMBERS.csv`: todas las características seleccionadas por MOFA/Otsu.
- `tables/30_LOGISTIC_FEATURE_IMPORTANCE_GLOBAL.csv`: **todas** las características de los 9 ajustes, con `importance_l2_standardized`, ranking global y dentro de vista, pesos MOFA, panel final 28 y coeficientes estandarizados de las 3 clases.
- `tables/31_LOGISTIC_FEATURE_IMPORTANCE_BY_CLASS.csv`: ranking de `|beta estandarizado|` separado para NP, FD y SO.
- `tables/32_MOFA_VS_LOGISTIC_RANK_CORRELATION.csv`: Spearman descriptivo (peso MOFA vs coeficiente logístico) por vista y modelo.
- `tables/33_PANEL_SIZES_AND_OVERLAPS.csv`: dimensión y solapamiento de paneles.
- `figures/IMPORTANCE_P{1,2,3}_{MOFA_OTSU,CORE18_47,FINAL28}_TOP30.png`.
- `figures/BETA_P{1,2,3}_{panel}_{NP,FD,SO}_TOP20.png`.
- `figures/ROC_OVR_...png` y `PR_OVR_...png` por conjunto y familia.

Los gráficos muestran un top 30/20 **solo para legibilidad**; no son umbrales de selección, las tablas conservan todos los predictores.

## Interpretación

Ridge multinomial con `alpha=0` y pesos inversos de frecuencia de clase. El coeficiente estandarizado se calcula como `beta (escala original) × SD(TRAIN)`. Importancia global = norma L2 de esos coeficientes entre las tres clases. Es un **ranking condicional al conjunto de predictores**, no una estimación causal, prueba de significación ni ranking estable garantizado frente a la colinealidad. Las importancias de modelos con paneles distintos no deberían interpretarse como magnitudes directamente comparables sin cautela.

Con TRAIN n=16 y FD n=3, las predicciones out-of-fold de TRAIN usan la `lambda` seleccionada en TRAIN completo, y la selección de biomarcadores ya depende de TRAIN; las métricas de TRAIN son optimistas. TEST n=10 **ya se ha usado previamente** y es un conjunto exploratorio reutilizado, no validación externa nueva. ROC/PR son one-vs-rest por clase y macro-promedios descriptivos.

**Versionado**: hacer `git add scripts/biomarker_selection/30_logistica_multinomial_paneles_y_pesos_MOFA.R scripts/biomarker_selection/30b_ROC_comparativa_modelos_OVR.R scripts/biomarker_selection/README_30.md` y no incluir eliminaciones locales de los scripts 28 y 29.

**Estado de comprobación**: se revisó estáticamente el código y se cotejaron los 6 modelos previos con el ZIP `20261031T223718_pid256178.zip`. No se ha ejecutado esta versión en R local desde el asistente.

## Nuevo v1.4.1 — comparacion ROC por clase SIN reentrenar

El script adicional `30b_ROC_comparativa_modelos_OVR.R` lee la ultima corrida
**completa** del script 30 (9 modelos: P1/P2/P3 × CORE18_47/FINAL28/MOFA_OTSU),
reutiliza los puntos de ROC y AUC del **mismo TEST** y dibuja las tres curvas
por clase y familia. Comprueba que el conjunto de muestras y etiquetas TEST
coincida entre los tres paneles. No vuelve a ajustar modelos ni a recalcular
probabilidades, lambda, biomarcadores o factores MOFA2.

Ejecutar directamente desde el repositorio, despues de una corrida v1.4.0 exitosa:

```bash
cd "$HOME/Documents/PhD/metatest/metatest"
set -o pipefail
Rscript --vanilla scripts/biomarker_selection/30b_ROC_comparativa_modelos_OVR.R \
  2>&1 | tee SCRIPT30_COMPARAR_ROC.log
printf 'EXIT_CODE=%s\n' "${PIPESTATUS[0]}"
```

Si tienes varias corridas y quieres escoger una concreta, usar antes del comando:
`MT2_30_ROC_SOURCE_DIR="/ruta/a/SCRIPT30_LOGISTICA_MOFA/fecha_pid"`.

Salidas **dentro de la corrida exitosa**:

- `tables/34_ROC_COMPARACION_TRES_MODELOS_OVR_TEST.csv`: 27 filas, tres modelos por familia y clase, con ROC-AUC, n positivos y n negativos.
- `figures/ROC_COMPARATIVAS_OVR/ROC_P1_NP_vs_rest_47_FINAL28_MOFA_TEST.png`, y analogas para P1/P2/P3 y NP/FD/SO (9 PNG).
- `figures/ROC_COMPARATIVAS_OVR/ROC_P1_NP_FD_SO_TRES_MODELOS_TEST.png` y analoga por familia (3 PNG).
- `figures/ROC_COMPARATIVAS_OVR/ROC_NP_P1_P2_P3_TRES_MODELOS_TEST.png` y analoga por clase (3 PNG).

Las curvas individuales previas siguen existiendo, sin modificaciones.
Al ejecutar de nuevo el script 30 **v1.4.1 completo**, esta comparacion se
genera automaticamente al terminar con los 9 estados PASS/PASS_REUSED.
**No necesitas reejecutar el script 30** para obtener las nuevas graficas
sobre tu corrida v1.4.0; basta con el `30b` independiente.

Cautela: TEST es pequeño y reutilizado; las diferencias de AUC no son pruebas
estadisticas de superioridad de los modelos.

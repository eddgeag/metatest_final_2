# METATEST2 — Script 30 (v1.1.0-repo)

**Ubicación versionada:** `scripts/biomarker_selection/30_logistica_multinomial_paneles_y_pesos_MOFA.R`.

**Ejecución:** análisis posterior a las etapas 18 y 28, sin reentrenamiento de MOFA2 ni repetición de PLS-DA, Elastic Net o sparse-group lasso. El análisis predictivo ajusta **regresiones multinomiales ridge nuevas**; no son modelos previamente entrenados. P1, P2, P3; paneles CORE28 (33), SHARED (exactamente 2 familias), SPECIFIC (1 familia), ALL, y CORE18_47 (47; consenso anterior). Un máximo de 15 ajustes; se marcan ausencias o fallos en `10_MODEL_STATUS.csv`.

El script busca automáticamente sus entradas dentro de la corrida local de CEDIA y registra las fuentes. **No ejecuta source() de los scripts anteriores**. Las salidas se guardan *fuera del repositorio Git*, bajo la ejecución de selección actual: `.../seleccion_mejor_escenario_train/<run>/SCRIPT30_LOGISTICA_MOFA/<fecha_pid>/`.

## Desde la raíz del repositorio Git

```bash
cd "$HOME/Documents/PhD/metatest/metatest"
python3 scripts/biomarker_selection/30_preflight.py
set -o pipefail
RUN="$HOME/Documents/PhD/metatest/results_cedia/metatest2_integral_66_20261027T013422Z"
mkdir -p "$RUN/script30_logs"
LOG="$RUN/script30_logs/script30_$(date +%Y%m%d_%H%M%S).log"
Rscript --vanilla scripts/biomarker_selection/30_logistica_multinomial_paneles_y_pesos_MOFA.R \
  2>&1 | tee "$LOG"
status=${PIPESTATUS[0]}
echo "R_EXIT_CODE=$status"
```

Para ajustar **solo CORE18_47**, sin repetir los 12 modelos:

```bash
MT2_30_ONLY_CORE18_47=1 Rscript --vanilla scripts/biomarker_selection/30_logistica_multinomial_paneles_y_pesos_MOFA.R
```

Se puede definir explícitamente `METATEST2_RUN_ROOT`, `METATEST2_SELECTION_RUN`, `METATEST2_SCRIPT28_EXPORT_DIR` y, cuando sea necesario, `METATEST2_RFM_P1`, `METATEST2_RFM_P2`, `METATEST2_RFM_P3`.

## Dependencias y procedencia

- Script 28: `FINAL_STORY_FIGURES_TABLES/<ultima>/00_GLOBAL/tables/FINAL_FEATURES_ALL.csv`, `FINAL_FEATURE_SHARING.csv`, `MOFA_NATIVE_IMPORTANCE.csv` y `EVIDENCE_A_ALL_MODEL_VOTES.csv`.
- Script 18: `BIOLOGICAL_PROGRAMS/<ultima>/FINAL_CONSENSUS_PANEL/<ultima>/tables/01_FINAL_PANEL_CORE_3OF3.csv` (debe contener exactamente 47 filas).
- Script 08/06: auditorías de programas y corroboración OBSERVED con la ruta a `ready_for_modeling.rds`, con vistas TX/PR/ME/CL y partición TRAIN/TEST preexistente.
- Script 05: `CROSSOMIC_BRIDGE/.../00_INPUT_AUDIT.csv` con ruta a MOFA2 RDS, **opcional** para exportar pesos firmados por factor.
- R: `glmnet` obligatorio; `MOFA2` opcional para pesos firmados. Para las normas L2 y el corte Otsu solo se usan los CSV de pesos exportados por el 28.

El script rellena `10_MODEL_STATUS.csv`, `11_MODEL_FEATURE_ELIGIBILITY.csv`, matrices y métricas de TRAIN-CV condicional y TEST exploratorio, coeficientes, modelos RDS, comparación CORE 33/47 y gráficos MOFA2. La transformación del peso MOFA es la norma L2 de factores (ya exportada); el umbral de Otsu sobre pesos al cuadrado es **descriptivo** y se calcula de forma idéntica por familia/vista. No es significación ni inferencia causal.

**Limitaciones:** TRAIN n=16, FD n=3, TEST reutilizado n=10. La selección de biomarcadores precede la CV; por tanto, no hay validación externa ni CV anidada de *todo* el pipeline. La comparación P2 con P3 no aísla FLI porque las selecciones de cuantiles difieren. El criterio FLI/PCR se rastrea en `04_FLI_PCR_TRACE.csv`.

**Para versionar**, añadir únicamente estos tres archivos al Git; NO añadir `results_cedia` ni el log. El archivo de salida `README.txt` registra la versión y el hash MD5 del código ejecutado junto con el commit Git observable al inicio. Guardar el commit después de haber actualizado el repositorio con el script 30, y después ejecutar el script.

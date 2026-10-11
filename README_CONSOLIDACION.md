# METATEST — consolidación local y ejecución (09-10-2026)

## Alcance y comprobaciones

Se compararon los cuatro scripts incluidos en `scripts/pipeline/` con la rama
`main` remota `eddgeag/metatest_final_2`: sus blobs Git son **idénticos**.
Por eso este paquete no los sobreescribe. Conserva el pipeline histórico,
las plantillas TSV, `run.sh` y todo el sistema de reproducibilidad previo.

Se incorporan los scripts de biomarcadores del ZIP sin modificar su contenido;
solo `SGL_FINAL_FIXED_V3.R` se **renombra** a
`16_sparse_group_lasso_hierarchical.R`: corresponde a la etapa 16.

**NO estaba en scripts.zip**: `scripts/selection/select_best_upstream_scenario_downstream.R`.
El instalador conserva el archivo que ya tienes en LOCAL. Comprueba que esté
presente antes de ejecutar la fase `selection`. Tampoco estaba en el ZIP
el directorio de datos ni los resultados de upstream; ambos se mantienen fuera
Git. Este paquete no demuestra que una ejecución completa haya terminado con éxito.

## Instalación SIN sobrescribir el histórico

1. Descomprime este ZIP **fuera** de la carpeta del repositorio.
2. Si ya estás en `analisis-nuevo`, preserva primero tu trabajo:
   `git status --short` y, opcionalmente, copia todo el repositorio como respaldo.
3. Ejecuta:

```bash
cd ~/Documents/PhD/metatest/metatest
bash /RUTA/EXTRAIDA/METATEST_CONSOLIDACION_20261009/instalar_en_repo.sh "$PWD"
git status --short
```

El instalador respalda los archivos locales sobrescritos bajo
`../metatest_consolidacion_backups/` y no borra ningún script de selección.

**Los dos scripts de selección marcados como eliminados en git status**:
`pca_espacios_train_project_test.R` y
`seleccionar_mejores_espacios_latentes.R` siguen disponibles en el historial
Git. Decide si sus eliminaciones son intencionales; no las restaures ni
confirmes automáticamente sin comprobarlo. Este paquete no las modifica.

## Orden de ejecución

```text
UPSTREAM (repositorio existente):
  scripts/run_upstream_only_local.sh
SELECCIÓN (script LOCAL no incluido en ZIP):
  scripts/selection/select_best_upstream_scenario_downstream.R
DESCUBRIMIENTO:
  01  bidirectional_framework_selected_scenarios
  02  latent_supervised_audit
  03  mofa_native_enrichment_part1
  04  class_pathway_bridge
  05  crossomic_bridge
  06  observed_train_corroboration
  07  test_transportability_audit
  08  define_biological_programs
  09  build_biological_synthesis_figures
  10  build_final_predictive_panel
  16  sparse_group_lasso_hierarchical (antes SGL_FINAL_FIXED_V3)
  17  evidenceA_plsda_pca_validation
  18  build_final_consensus_panel
  19  probe_final_small_panel_console (opcional, solo terminal)
POSTANÁLISIS, RStudio, NO automáticos:
  19_probe_final_small_panel_SCREEN_PLOTS.R (interactivo; pide ENTER)
  20_probe_small_final_panel_RSTUDIO.R
  21_close_project_final_panel_summary_RSTUDIO.R
  22_final_panels_by_question_RSTUDIO.R
```

`18` necesita la salida de `17` y puede usar la de `16` como anotación.
La ejecución secuencial deja que cada script encuentre el `LATEST_RUN.txt`
del anterior. El runner congela `SELECTION_RUN` durante la misma ejecución.

## Ejecutar sobre tus resultados existentes (local)

```bash
cd ~/Documents/PhD/metatest/metatest
ROOT="$HOME/Documents/PhD/metatest/results_upstream_local_20260908"
bash scripts/orchestration/run_metatest_vigente.sh \
  --phase biomarkers --results-root "$ROOT" --dry-run

# Cuando el entorno esté preparado y renv.lock actualizado:
bash scripts/orchestration/run_metatest_vigente.sh \
  --phase biomarkers --results-root "$ROOT" \
  2>&1 | tee "$ROOT/orquestacion_vigente.log"
```

Si los resultados previos ya contienen etapas completas, reanuda con
`--from 05` (o la etapa pertinente). No uses `--from` si las etapas anteriores
no están completas. El runner NO reutiliza ni descarta outputs automáticamente:
los propios scripts gestionan sus salidas. Para acotar una ejecución,
`--from 01 --to 09`, por ejemplo.

## Ejecución desde cero en un servidor

```bash
git clone https://github.com/eddgeag/metatest_final_2.git metatest
cd metatest
# Preparar ../data con el mismo dataset y sus hashes de entrada.
# Preparar R 4.4.2, Renv, CmdStan 2.36.0, Basilisk/MOFA2.
# Ver instrucciones de entorno más abajo.

bash scripts/orchestration/run_metatest_vigente.sh --phase all --dry-run
bash scripts/orchestration/run_metatest_vigente.sh --phase all
```

`--phase all` recorre el upstream COMPLETO de las plantillas: puede incluir
201 escenarios. `--phase selection` solo corre el selector existente en LOCAL.
Para no repetir barridos pesados, usa `--phase biomarkers` con el
`--selection-run /RUTA/EXACTA/` de la corrida actual.

Si copias RESULTADOS de otra máquina, sus CSV y `LATEST_RUN.txt` pueden contener
rutas absolutas del equipo anterior. El runner y los scripts NO reescriben esas
referencias. Recalcula la selección en la nueva máquina o verifica y remapea
las rutas de entradas y manifests antes de iniciar. `git clone` NO descarga
resultados, datos ni CmdStan, y `git pull` no arregla rutas antiguas.

## Entorno R y bloqueo de paquetes (REQUISITO ANTES DE PUBLICAR)

La rama `main` contiene `renv.lock` con R 4.4.2 / Bioconductor 3.20,
pero **no contiene** `msigdbr`, `msgl` ni `sglOptim`. No es correcto decir que
el nuevo pipeline se reproduce solo con el `renv.lock` histórico.
En el portátil que ya ejecutó el pipeline:

```bash
cd ~/Documents/PhD/metatest/metatest
Rscript scripts/reproducibility/actualizar_renv_vigente.R
# Verificar y revisar siempre qué cambió:
git diff -- renv.lock
Rscript -e 'source("renv/activate.R"); print(renv::status())'
```

`actualizar_renv_vigente.R` restaura primero el lock existente e instala las
versiones históricas `sglOptim 1.3.8` y `msgl 2.3.9`, además de
`msigdbr 26.1.1`. Sus paquetes deben poder instalarse desde los archivos
de CRAN: si un mirror falla, revisar los logs, no modificar las versiones
silenciosamente. La actualización realiza `renv::snapshot(type="all")`:
**revisar el diff** y versionar `renv.lock` solo tras comprobar los cambios.
No se modifica el `renv.lock.validated` del release de agosto.

Para la ejecución activa el entorno con
`R_PROFILE_USER=scripts/reproducibility/activate_renv.Rprofile` y
`METATEST_REPO_ROOT=$(pwd)`, ya configurados por el orquestador.
El script de preflight detiene la ejecución si el lock y los paquetes instalados
no concuerdan. Configurar CmdStan 2.36.0 con `CMDSTAN=/ruta/.../cmdstan-2.36.0`
si no está en `~/.cmdstan/` y la capa Basilisk/MOFA2, por separado.

**MSigDB** se obtiene mediante `msigdbr` y puede requerir descarga externa;
el repositorio no incluye estos datos. Una prueba de conexión y disponibilidad
de gene sets es necesaria incluso después de instalar el paquete R.

## Docker opcional (después de actualizar renv.lock)

El Dockerfile es una base de construcción, NO una imagen validada en este chat.
Antes de construir, habilitar el seguimiento de `Dockerfile.metatest-vigente`
y `.dockerignore`, y congelar el nuevo `renv.lock`. Por defecto `docker build`
usa `.dockerignore`; el instalador copia `.dockerignore.metatest-vigente` como
`.dockerignore` al aplicar el paquete. El container tiene que compilar
paquetes R, CmdStan y preparar el backend Python de MOFA2; la primera
construcción necesita Internet y puede ser costosa.

```bash
docker build -f Dockerfile.metatest-vigente -t metatest-vigente:20261009 .
mkdir -p ../results
docker run --rm -it --cpus=4 --memory=32g \
  -v "$(realpath ../data):/workspace/data:ro" \
  -v "$(realpath ../results):/workspace/results" \
  metatest-vigente:20261009 \
  bash scripts/orchestration/run_metatest_vigente.sh --phase biomarkers
```

El ejemplo asume que `../results` contiene resultados upstream y selección
válidos **creados o remapeados en esas rutas**. No use la imagen como evidencia
de reproducibilidad hasta verificar instalación, CmdStan, Basilisk, inputs,
y una corrida mínima. Para mayor reproducibilidad fijar el digest SHA256 de la
imagen base y de la imagen final; conservar manifiestos de datos y entorno.

## Git: consolidación inicial y posteriores cambios

Se recomienda incorporar el trabajo mediante la rama `analisis-nuevo`, mantener
`main` como rama de ejecución y conservar el tag histórico.

```bash
cd ~/Documents/PhD/metatest/metatest
git fetch origin
git switch analisis-nuevo
git add -A -- scripts/biomarker_selection scripts/orchestration scripts/reproducibility scripts/plots_presentation
git add -- scripts/selection/select_best_upstream_scenario_downstream.R renv.lock Dockerfile.metatest-vigente .dockerignore .gitignore README_CONSOLIDACION.md
# Comprueba el contenido exacto ANTES de confirmar:
git diff --cached --stat
git diff --cached --name-status
git commit -m "Integrar pipeline bidireccional y descubrimiento reproducible"
git merge origin/main
# Verificar que los scripts están y los tests de entorno funcionan.
git switch main
git merge --ff-only analisis-nuevo
git push origin main
```

Si `main` no es ancestro de `analisis-nuevo` o hay conflictos, resolver de forma
explícita; no forzar el push. Las eliminaciones de los dos scripts antiguos de
selección NO se añaden con el `git add` anterior. Si son intencionales,
incluirlas por separado con `git add -u -- scripts/selection/`.

Para actualizaciones futuras: `git pull` DESCARGA, `git push` PUBLICA
commits. **No crean commits de modificaciones locales**. Usar:

```bash
git pull --ff-only
# editar scripts de figuras en scripts/plots_presentation/
bash scripts/orchestration/publicar_codigo.sh "Mejorar figuras de tesis"
```

El publicador prepara solamente scripts y docs en carpetas permitidas, crea
el commit, hace `git pull --rebase origin main` y publica. No agrega datos ni
resultados. Alternativamente, usar git add/commit/pull/push manualmente.

## Límites de la auditoría realizada

Verificados: integridad de los 4 scripts `pipeline`, orden 16, lectura estática
de dependencias, nombres y rutas, sintaxis Bash del orquestador.
NO verificados: ejecución R (este contenedor de preparación no dispone de R),
instalación completa de renv/CmdStan/Basilisk ni coherencia de los artefactos
locales reales (no se han suministrado). El selector LOCAL no estaba en el ZIP.

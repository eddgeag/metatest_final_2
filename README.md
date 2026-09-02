# METATEST

Pipeline de comparacion de modelamiento2 y modelamiento3. La rama de trabajo
es `main`; el tag `metatest-reproducible-20260830` conserva el release historico.
No contiene datos, resultados, paquetes R ni imagenes de ejecucion.

## Organizacion

Un directorio de trabajo con `metatest/` (este repositorio), `data/` (inputs),
`assets/` (runtime y backend sellados) y `results/` (salidas).
Las rutas explicitas son relativas al directorio donde se invoca `run.sh`.
Los valores por defecto son los directorios hermanos del repositorio.

## Planificar sin ejecutar

```bash
bash run.sh --list
bash run.sh --backend enroot --stan-o 1 --data ../data --assets ../assets \
  --outdir ../results/barrido_completo_O1 --check
```

El plan completo contiene 201 escenarios: 67 por plantilla, incluido un q00
por familia. Cada escenario ejecuta ambos modelos. `--list` y `--dry-run`
no crean salidas. `--check` verifica el checkout y las rutas sin ejecutar R ni
enviar trabajos; no sustituye al smoke y los hashes completos dentro del job.

## Ejecutar cuando se autorice el barrido

```bash
bash run.sh --backend enroot --stan-o 1 --data ../data --assets ../assets \
  --outdir ../results/barrido_completo_O1 --slurm
```

`run.sh` usa `run.sbatch`: una tarea por escenario y una tarea simultanea
(`%1`). Valores por defecto: cpu-dev, 8 CPU, 64G y 8 horas por tarea.
Son solicitudes de recursos, no garantias de duracion ni de cuota disponible.
Al enviar se imprimen el ID, las rutas de logs y el comando de recoleccion.
No ejecutar calculos en login. Puede cerrarse la conexion despues del envio.

Un solo escenario, incluido el base:

```bash
bash run.sh --backend enroot --stan-o 1 \
  --scenario P1_ALL_qtx00_qpr00_qme00_qcl00 \
  --data ../data --assets ../assets --outdir ../results/prueba_q00 --slurm
```

Una fila (1=q00; no cuenta comentarios/cabecera):

```bash
bash run.sh --backend enroot --stan-o 1 \
  --plantillas templates/current/plantilla_2_con_FLI_barrido.tsv --fila 2 \
  --data ../data --assets ../assets --outdir ../results/prueba_fila --slurm
```

`--plantillas` y `--scenario` admiten varias apariciones. Sin `--slurm`, el
backend Enroot se ejecuta en el nodo actual: usar solo un nodo de computo.
`--backend local` selecciona y comprueba el entorno nativo descrito abajo;
no instala ni certifica equivalencia cientifica con Enroot. El rebuild automatico
O1 solo esta disponible en Enroot. O0 sigue
siendo el valor por defecto; O1 debe pedirse explicitamente.

## Portatil o servidor propio sin Enroot, Docker ni Slurm

Se requieren R 4.4.2, los 302 paquetes y versiones de `renv.lock`, CmdStan
2.36.0, y el entorno Basilisk de MOFA2 1.16.0 con mofapy2 0.7.0 ya instalado.
Git no transporta estas dependencias. Un clon nuevo necesita restaurar su
biblioteca; no basta `git clone` para preparar una maquina vacia.
No usar el antiguo restore de CEDIA como instalador nativo automatico.

El runner busca una unica biblioteca de proyecto `renv/library/*/R-4.4/...`.
Opcionalmente `METATEST_R_LIBRARY` permite elegir explicitamente otra
biblioteca ya preparada y verificable. `METATEST_CMDSTAN` (o `CMDSTAN`)
selecciona CmdStan; el valor por defecto es `~/.cmdstan/cmdstan-2.36.0`.
Basilisk consulta su propia cache, o `BASILISK_EXTERNAL_DIR` si se especifica.
No hay rutas de usuario incrustadas en el codigo. No se descarga Python ni
se restauran paquetes durante una corrida: si faltan, la comprobacion se detiene.

Rscript y sus workers reciben la biblioteca del proyecto mas la biblioteca
base de R. Se ignoran los perfiles y bibliotecas globales, no se modifica
`~/.Rprofile`, y se fija `RETICULATE_PYTHON` al Python del entorno MOFA.
Se conserva la cache de Basilisk (incluidos otros entornos); sus locks necesitan
permisos de escritura. El CmdStan instalado no se limpia ni se reconstruye
automaticamente. La compilacion normal de modelos puede generar su cache.

Primero, solo dependencias e inputs (sin pipeline y sin compilar Stan):

```bash
bash scripts/portable/check_native.sh --data ../data --outdir ../results/prueba_nativa
```

El smoke verifica versiones/RemoteSha registrados, carga paquetes incluido
igraph y prueba la conexion R-reticulate-Basilisk-MOFA. Guarda logs y una huella
del runtime fuera del repositorio. No valida rendimiento, toda la compilacion
C++, ni resultados cientificos. Repite estas comprobaciones al ejecutar.

Despues de `NATIVE_CHECK=PASS`, un unico escenario q90 con ambos modelos:

```bash
bash run.sh --backend local --stan-o 0 \
  --scenario P1_TODAS_qtx90_qpr45_qme90_qcl45 \
  --data ../data --outdir ../results/prueba_nativa_P1_q90
```

No anadir `--slurm` en una maquina sin Slurm. La corrida es en primer plano;
no cerrar la terminal durante la prueba. Se usan los limites conservadores
del runner (2 cores de pipeline, 1 interno/BLAS), no todos los CPUs detectados.
El escenario q90 no sustituye una aceptacion q00 nativa posterior.
Sin `--scenario`, `--fila` ni `--pilot`, `run.sh` ejecuta los 201 escenarios.

La huella nativa incluye versiones/builds de paquetes, binarios criticos,
inventario Python, rutas y compilador. `--resume` rechaza una huella distinta
o ausente; se debe usar un outdir nuevo si cambia el entorno. No es una imagen
sellada ni garantiza identidad binaria de todo el sistema operativo.

## Reproducibilidad y trabajo diario

- Referencias: job 24178 (q00 reanudado valido), 24294 (q00 limpio aceptado).
- Imagen/backend: identidades y hashes en `config/cedia_release.env`.
  Estos archivos pesados se suministran fuera de Git. No basta clonar para
  obtenerlos ni para obtener los datos privados.
- Enroot usa `create` y `start`, sin squashfuse ni Rscript en el host.
  Los datos y el backend se montan RO; el estado Basilisk es efimero RW.
- O1 reconstruye solo el CmdStan del contenedor temporal, no la imagen sellada.
  Los pilotos pasaron y sus 36 CSV de resumen q90 coincidieron con O0; esto no
  demuestra equivalencia de todos los artefactos ni de los 201 escenarios.
- Cada corrida registra commit de codigo, base del runtime, manifiesto,
  escenario, perfil O0/O1, tiempos y validaciones. El barrido O1 sigue siendo
  una corrida portable experimental, no un nuevo release cientifico sellado.
- `main` puede recibir commits de analisis post hoc (`scripts/selection/` u
  otro directorio). El runner exige checkout limpio, descendencia del release
  y nucleo cientifico/runtime sin cambios respecto a ese release. Cambiar el
  nucleo requiere una nueva validacion; no eliminar ese control para avanzar.
- No modificar el checkout mientras haya tareas pendientes: cada tarea
  comprueba la instantanea de codigo tomada al enviar el array.
- Usar un outdir nuevo para cada corrida limpia. `--resume` solo acepta el mismo
  perfil; no usarlo para convertir un benchmark previo en una corrida limpia.
- El historial y los runners bajo `scripts/reproducibility/` documentan el
  release original. La entrada de trabajo actual es `run.sh`; no usar el antiguo
  `LANZAR_BARRIDO_COMPLETO_O1_20260901.sh` tras consolidar o cambiar las rutas.

```bash
git pull --ff-only
# Editar, revisar y verificar el codigo post hoc.
git add scripts/selection/archivo.R
git commit -m "Describe el cambio post hoc"
git push origin main
```

No usar `git add .` indiscriminadamente. `renv/library/` es una dependencia
local ignorada, no basura. `reference_maps/` es generado; el orquestador lo
reconstruye desde los scripts de `scripts/utilities/`.

## Recoger resultados

```bash
bash scripts/portable/collect.sh ../results/barrido_completo_O1
```

El recolector solo permite recopilar cuando no hay tareas activas del array y empaqueta logs y
tablas de resumen, no inputs ni RDS. Un PASS funcional no elimina advertencias
estadisticas: por ejemplo, P3 q90/m2 conserva un Pareto k maximo cercano a 0.814
ya presente con O0, que requiere revisar PSIS-LOO antes de interpretarlo.

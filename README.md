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
`--backend local` usa el R/CmdStan ya configurado; no instala ni certifica ese
entorno. El rebuild automatico O1 solo esta disponible en Enroot. O0 sigue
siendo el valor por defecto; O1 debe pedirse explicitamente.

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

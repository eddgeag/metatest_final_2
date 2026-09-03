# METATEST

Pipeline reproducible de comparación de `modelamiento2` y `modelamiento3`.
La rama de trabajo es `main`; el tag `metatest-reproducible-20260830` conserva
el release histórico. La entrada común es `run.sh`, tanto para ejecución
nativa como para Enroot. `run.sbatch` es su lanzador para Slurm.

El barrido completo contiene **201 escenarios: 67 por plantilla, incluido
un q00 base en cada una**. Cada escenario ejecuta ambos modelamientos.
Los ejemplos de pruebas y de barrido son alternativas: no hay que lanzarlos
todos. Una corrida limpia siempre utiliza una carpeta de salida nueva.

## 1. Repositorio, datos y entorno

Todos los comandos de este documento se ejecutan desde la raíz del repositorio.
Para un clon nuevo, usar un clon completo con el historial y el tag histórico;
el control de procedencia necesita ambos.

```bash
git clone --branch main https://github.com/eddgeag/metatest_final_2.git metatest
cd metatest
```

Organización recomendada:

| Carpeta | Contenido |
| --- | --- |
| `metatest/` | Este repositorio y sus scripts versionados. |
| `data/` | Entradas privadas; carpeta hermana del repositorio. |
| `results/` | Resultados, logs y evidencias; carpeta hermana del repositorio. |
| `assets/` | Imagen y backend sellados; solo se necesitan con Enroot. |

`--data`, `--outdir` y `--assets` aceptan rutas absolutas o relativas al
directorio desde el que se invoca `run.sh`. Si se omiten, se usan las carpetas
hermanas del repositorio. `--outdir` debe quedar fuera del repositorio, los
datos y los assets. Se recomienda evitar espacios en las rutas del proyecto
y de CmdStan, por las restricciones de Make.

Git **no incluye** los datos privados, los paquetes R, CmdStan, los entornos
Python ni las imágenes. Clonar no prepara automáticamente una máquina vacía.
Las entradas deben coincidir con `reproducibilidad/release/inputs.sha256`;
cambiar su ubicación no equivale a autorizar otro dataset.

| Backend | Requisitos antes de ejecutar |
| --- | --- |
| `local` | Linux/Bash 4+, R 4.4.2, biblioteca con los 302 paquetes de `renv.lock`, CmdStan 2.36.0, entorno Basilisk de MOFA2 1.16.0 con mofapy2 0.7.0, Python, compilador C++, Make y utilidades del sistema. No necesita Docker, Enroot ni Slurm. |
| `enroot` | Host compatible con el runtime sellado, Enroot, `unsquashfs`, utilidades del sistema y los dos assets sellados. No requiere `Rscript` en el host ni `squashfuse`: utiliza `enroot create` y después `enroot start`. |
| Cualquiera con `--slurm` | Slurm y acceso permitido a partición, cuenta y recursos. Los requisitos del backend deben estar disponibles en los nodos de cómputo. |

Los dos archivos de `assets/` son:

- `metatest-runtime-metatest-reproducible-20260830.sqsh`
- `metatest-mofa-backend-metatest-reproducible-20260830.tar`

Sus identidades están registradas en `config/` y `reproducibilidad/release/`.
Enroot verifica los hashes al ejecutar; monta datos y payload del backend en
solo lectura y mantiene el estado de Basilisk en una ubicación efímera escribible.

## 2. Plantillas y selección de escenarios

| Familia | Archivo en `templates/current/` | Escenarios |
| --- | --- | ---: |
| P1: todas las clínicas | `plantilla_1_todas_clinicas_barrido.tsv` | 67, con q00 |
| P2: con FLI | `plantilla_2_con_FLI_barrido.tsv` | 67, con q00 |
| P3: sin FLI | `plantilla_3_sin_FLI_barrido.tsv` | 67, con q00 |

Consultar planes sin ejecutar ni crear salidas:

```bash
# Las tres plantillas, con sus tres escenarios base: 201 filas de análisis.
bash run.sh --list

# Los tres q90 balanceados: un escenario por plantilla, sin q00.
bash run.sh --pilot q90 --list

# Una sola plantilla completa: 67 escenarios.
bash run.sh --plantillas templates/current/plantilla_2_con_FLI_barrido.tsv --list

# Primera fila de análisis de P1: q00; no cuenta comentarios ni cabecera.
bash run.sh --plantillas templates/current/plantilla_1_todas_clinicas_barrido.tsv \
  --fila 1 --list
```

`--plantillas` y `--scenario` se pueden repetir. `--fila N` exige una sola
plantilla y no se combina con `--scenario` ni `--pilot`. No editar las
plantillas para seleccionar filas: utilizar estas opciones.

Sin `--plantillas`, se seleccionan las tres plantillas actuales. Sin
`--scenario`, `--fila` ni `--pilot`, se ejecutan todas las filas seleccionadas,
incluidos los q00. `--dry-run` muestra las rutas resueltas y el plan sin ejecutar.

## 3. Ejecutar en local, sin contenedor ni Slurm

La biblioteca nativa y los demás requisitos deben estar preparados previamente.
No se instalan ni actualizan dependencias durante una corrida. Se busca una
única biblioteca del proyecto bajo `renv/library/*/R-4.4/`; opcionalmente,
`METATEST_R_LIBRARY` permite elegir otra biblioteca ya preparada.
`METATEST_CMDSTAN` o `CMDSTAN` seleccionan CmdStan; el valor por defecto es
`~/.cmdstan/cmdstan-2.36.0`. Basilisk utiliza su caché o `BASILISK_EXTERNAL_DIR`.

El runner aísla las bibliotecas y perfiles globales de R, configura el Python
de MOFA y comprueba una huella del entorno. No modifica `~/.Rprofile` ni
reconstruye automáticamente el CmdStan del host. La caché de Basilisk necesita
permisos de escritura para sus locks. No copiar bibliotecas compiladas entre
sistemas distintos suponiendo que son intercambiables.

Primero, con el checkout sin cambios pendientes:

```bash
sha256sum -c PORTABLE_SHA256SUMS.txt

bash run.sh --backend local --stan-o 0 \
  --data ../data --outdir ../results/barrido_local_O0 --check

bash scripts/portable/check_native.sh \
  --data ../data --outdir ../results/comprobacion_nativa
```

`--check` es estático: no ejecuta R, no comprueba todo el runtime ni envía jobs.
`check_native.sh` verifica inputs, dependencias y el puente
R-reticulate-Basilisk-MOFA, pero **no ejecuta el pipeline ni prueba la
compilación de los modelos Stan**. Versiones coincidentes no garantizan por
sí solas que todos los paquetes fueron compilados de forma compatible.

Después de `NATIVE_CHECK=PASS`, prueba limpia de un escenario con ambos modelos:

```bash
bash run.sh --backend local --stan-o 0 \
  --scenario P1_TODAS_qtx90_qpr45_qme90_qcl45 \
  --data ../data --outdir ../results/prueba_local_P1_q90_O0
```

Barrido completo, cuando se decida ejecutarlo:

```bash
bash run.sh --backend local --stan-o 0 \
  --plantillas templates/current/plantilla_1_todas_clinicas_barrido.tsv \
  --plantillas templates/current/plantilla_2_con_FLI_barrido.tsv \
  --plantillas templates/current/plantilla_3_sin_FLI_barrido.tsv \
  --data ../data --outdir ../results/barrido_local_O0
```

Omitir las tres opciones `--plantillas` produce el mismo plan. `./run.sh`
también es válido si se conserva su permiso de ejecución. La corrida nativa
es en primer plano; no cerrar la terminal ni suspender el equipo mientras
esté ejecutándose. No añadir `--slurm` en una máquina sin Slurm.

El perfil actual limita el pipeline a 2 cores, el paralelismo interno de
BRMS/BLAS a 1 y Make a `-j1`. Solicitar más CPUs no cambia automáticamente
esos límites ni convierte el barrido en ejecución paralela.

## 4. Ejecutar en un servidor con Slurm

### Consultar los límites vigentes

Antes del envío, consultar el propio scheduler; los nombres de las particiones
y sus límites no son universales:

```bash
sinfo -o '%P|estado=%a|max_time=%l|nodos=%D|cpu_nodo=%c|memoria_MB_nodo=%m'
scontrol show partition
sacctmgr show assoc where user="$USER" \
  format=Cluster,Account,User,Partition,QOS,GrpJobs,MaxJobs,MaxSubmit,GrpTRES,MaxTRES,GrpWall,MaxWall -P
sacctmgr show qos -P
scontrol show config | awk '/MaxArraySize|MaxJobCount|AccountingStorageEnforce/'
```

Comprobar `MaxTime`, `MaxWall`, CPUs/memoria, trabajos simultáneos, máximo de
trabajos enviados y tamaño permitido del array. La cuenta, su jerarquía y las
QoS pueden imponer restricciones adicionales. Si no se permite consultarlas,
confirmarlas con la administración; una consulta vacía no acredita ausencia
de límites. Véase la [jerarquía de límites de Slurm](https://slurm.schedmd.com/resource_limits.html).

El runner incorpora como valores por defecto `cpu-dev`, 8 CPU, 64G y 8 horas.
La partición es un valor heredado: elegir explícitamente una admitida por el
servidor. Definir estas dos variables en la misma terminal que hará el envío:

```bash
METATEST_PARTITION='NOMBRE_PARTICION'
METATEST_WALLTIME='08:00:00'
```

Sustituir `NOMBRE_PARTICION` por el nombre real. Elegir el tiempo a partir del
escenario más lento medido, con margen para preparación, compilación y cierre.
Las 8 horas son una solicitud inicial, **no una garantía para cada escenario**.

| Valor de `METATEST_WALLTIME` | Límite solicitado por escenario |
| --- | --- |
| `08:00:00` | 8 horas |
| `1-00:00:00` | 24 horas, si partición y cuenta/QoS lo permiten |
| `5-00:00:00` | 5 días, solo si los límites efectivos permiten esa duración |

Si la partición informa `MaxTime=5-00:00:00`, los cinco días son un techo de
esa partición, no una duración recomendada ni el límite acumulado de todo el
barrido. Pedir más tiempo no acelera el programa. Al agotarlo, Slurm termina
el job; aumentar el valor requiere otra asignación autorizada.
Véase [`sbatch --time`](https://slurm.schedmd.com/sbatch.html#OPT_time).

### Comprobación y envío de las tres plantillas

Este ejemplo utiliza Enroot y una carpeta nueva para el barrido completo:

```bash
bash run.sh --backend enroot --stan-o 0 \
  --data ../data --assets ../assets \
  --outdir ../results/barrido_enroot_O0 --check
```

Solo después de revisar el plan, los recursos y el `CHECK=PASS`, enviar:

```bash
bash run.sh --backend enroot --stan-o 0 \
  --plantillas templates/current/plantilla_1_todas_clinicas_barrido.tsv \
  --plantillas templates/current/plantilla_2_con_FLI_barrido.tsv \
  --plantillas templates/current/plantilla_3_sin_FLI_barrido.tsv \
  --data ../data --assets ../assets \
  --outdir ../results/barrido_enroot_O0 \
  --partition "${METATEST_PARTITION:?Define una particion admitida}" \
  --cpus 8 --mem 64G --time "${METATEST_WALLTIME:?Define el limite por escenario}" \
  --slurm
```

El envío utiliza **el `run.sbatch` ya versionado**: array `0-200%1`, un
escenario con ambos modelos por tarea y como máximo una tarea activa a la vez.
No crea 201 tareas simultáneas. Cada tarea recibe su propio límite de tiempo;
el calendario total puede superar cinco días aunque ninguna tarea los supere.
`%1` no garantiza el orden de los índices ni evita los límites de trabajos
pendientes/enviados. Véase [arrays de Slurm](https://slurm.schedmd.com/job_array.html).

Si la cuenta no admite 201 tareas enviadas, seleccionar una plantilla por
envío (67 tareas) o un subconjunto menor con `--scenario`, respetando su cuota.
No enviar varios arrays a la misma salida simultáneamente. El límite de una
tarea activa es por array; no aumenta una cuota de usuario de un solo job.

Para probar únicamente los tres q90 antes del barrido, usar `--pilot q90`
con otro `--outdir`; para un único q00, usar
`--scenario P1_ALL_qtx00_qpr00_qme00_qcl00`. Las opciones de selección de la
sección 2 se aplican también a Slurm.

Slurm también puede ejecutar el backend nativo: utilizar `--backend local
--stan-o 0`, omitir `--assets` y escoger otra salida. La biblioteca R,
CmdStan, Python y los datos deben estar accesibles en el nodo asignado. No
ejecutar la comprobación nativa pesada ni los modelos en un nodo de login.

Para enviar directamente **un solo escenario** con el batch, como alternativa
al array gestionado por `run.sh`:

```bash
mkdir -p ../results/prueba_batch_q00_O0/run_reports
sbatch --partition "${METATEST_PARTITION:?}" \
  --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=64G \
  --time "${METATEST_WALLTIME:?}" --no-requeue \
  --output=../results/prueba_batch_q00_O0/run_reports/slurm_%j.out \
  --error=../results/prueba_batch_q00_O0/run_reports/slurm_%j.err \
  run.sbatch --backend enroot --stan-o 0 \
  --scenario P1_ALL_qtx00_qpr00_qme00_qcl00 \
  --data ../data --assets ../assets --outdir ../results/prueba_batch_q00_O0
```

Hacer este envío desde la raíz del repositorio. No añadir `--slurm` después
de `run.sbatch`, no ejecutar `bash run.sbatch` en login y no enviar todo el
barrido como un único job monolítico. El envío directo no guarda la misma
instantánea de código que el array gestionado; para el barrido usar
`bash run.sh ... --slurm`.

Una vez aceptado el envío se puede cerrar SSH. No modificar ni actualizar el
checkout hasta que hayan terminado las tareas pendientes y en ejecución;
estas utilizan los archivos de ese checkout.

## 5. O0, O1 y tiempo de ejecución

`--stan-o 0` y `--stan-o 1` seleccionan el nivel `O` de compilación **C++**
de CmdStan. No son el número de CPUs, cadenas o iteraciones, ni la opción
`--O1` del traductor `stanc`.

| Perfil | Disponibilidad y efecto en este runner |
| --- | --- |
| O0 (`--stan-o 0`) | Valor por defecto; referencia común para nativo y Enroot. |
| O1 (`--stan-o 1`) | Solo Enroot, experimental. Reconstruye CmdStan dentro del contenedor temporal de cada escenario; no cambia los assets sellados ni el CmdStan del host. |

O1 puede reducir el tiempo de muestreo, pero añade trabajo de compilación y
no asegura menor tiempo total. La [documentación de GCC](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html)
explica este intercambio; no es una promesa de aceleración para este pipeline.
El runner no admite `--stan-o 2` ni `--stan-o 3`. Reservar más CPU o más días
tampoco cambia por sí solo el paralelismo interno.

Para evaluar O1, comparar escenarios iguales, desde cero, con los mismos
datos, commit, hardware y recursos, usando salidas separadas O0/O1. Incluir
preparación y reconstrucción de CmdStan en el tiempo total, además del
muestreo. Revisar diagnósticos y resultados: un PASS funcional no certifica
equivalencia estadística ni identidad de todos los artefactos. Mantener O0
para la comparación directa con el perfil nativo actual.

## 6. Logs, finalización y reanudación

El envío mediante `run.sh --slurm` imprime el ID del array, los comandos de
consulta y `tail -F`, y la ruta de recolección. Cada envío deja sus logs bajo
`OUTDIR/run_reports/submission_.../slurm_<array>_<tarea>.out` y `.err`.
Cada escenario guarda petición, procedencia, tiempos por fase y validaciones
en su propio informe dentro de `OUTDIR/run_reports/`.

```bash
squeue -u "$USER"
```

Cuando no quede ninguna tarea activa, recoger la evidencia del barrido:

```bash
bash scripts/portable/collect.sh ../results/barrido_enroot_O0
```

Para local, sustituir la ruta por el `--outdir` utilizado. El recolector
empaqueta logs y tablas de resumen, no inputs ni RDS. Una cola vacía o un
`SUBMITTED=PASS` no acreditan éxito: comprobar el estado final de cada tarea
en `sacct`, los informes de validación y que estén todos los escenarios
seleccionados. Una tarea fallida no cancela automáticamente las demás del array.

Un `TIMEOUT` puede dejar pasos completos aprovechables, pero no significa
que el escenario terminó. No hay reenvío automático ni continuación exacta
de una cadena MCMC interrumpida. Para reanudar, revisar primero los informes
y añadir `--resume` al mismo comando, conservando código, datos, rutas,
backend, O0/O1 y, en nativo, la huella del runtime. Se pueden seleccionar
solo los escenarios fallidos. No modificar marcas ni huellas para forzarlo.
Una reanudación no es una medición desde cero; para una corrida limpia, usar
otro `--outdir` y omitir `--resume`.

## 7. Reproducibilidad y trabajo normal con Git

- Referencias históricas: job 24178, q00 reanudado válido; job 24294, q00
  limpio aceptado. Son identificadores de evidencia, no requisitos para
  volver a ejecutar ni versiones del código.
- El runner exige cambios tracked confirmados, manifiesto íntegro,
  descendencia del release y núcleo científico/runtime sin alteraciones
  respecto al release. Cambiar ese núcleo requiere nueva validación.
- Los análisis post hoc pueden añadirse, por ejemplo, a `scripts/selection/`.
  No es necesario crear otra rama para el flujo de trabajo habitual en `main`.
- La huella nativa registra versiones, builds, binarios relevantes, rutas y
  compilador. Compartir commit y dataset no equivale a compartir todos los
  binarios del sistema. No se exige la misma huella entre hosts distintos.
- Mantener las advertencias estadísticas separadas del éxito técnico: revisar
  convergencia, PSIS-LOO/Pareto k y limitaciones de tamaño muestral antes de
  interpretar modelos. Un piloto q90 no valida los 201 escenarios ni sustituye q00.
- `renv/library/` es una dependencia local ignorada, no basura. Los mapas de
  referencia se regeneran desde `scripts/utilities/`. No añadir datos,
  resultados, cachés o assets pesados al repositorio.

Antes de editar, sin corridas pendientes o activas que usen el checkout:

```bash
git fetch origin
git pull --ff-only
git status --short --branch
```

Después de editar únicamente este documento:

```bash
git diff --check
git diff -- README.md
git add -- README.md
git commit -m "Documenta ejecucion local y Slurm de las tres plantillas"
git push origin main
```

`README.md` no forma parte de `PORTABLE_SHA256SUMS.txt`; una modificación
exclusivamente documental no exige regenerar ese manifiesto. Sí crea un
commit nuevo: no mezclarlo con una reanudación cuyo perfil exige el commit
anterior. No usar `git add .` indiscriminadamente ni forzar un push.

En el otro clon, con el checkout limpio y sin corridas activas:

```bash
git fetch origin refs/heads/main:refs/remotes/origin/main
git pull --ff-only origin main
git symbolic-ref HEAD
git rev-parse HEAD
git rev-parse origin/main
git rev-list --left-right --count HEAD...origin/main
git status --short --branch
sha256sum -c PORTABLE_SHA256SUMS.txt
```

Ambos clones deben mostrar `refs/heads/main`, el mismo `HEAD`, comparación
`0 0`, ningún cambio pendiente y todos los hashes correctos. El `fetch` con
destino explícito actualiza también `origin/main` en clientes antiguos que
pueden dejarlo desactualizado tras un `pull origin main`; no reescribe el
historial. Véase [refspecs de Git](https://git-scm.com/book/en/v2/Git-Internals-The-Refspec).

Tras sincronizar, repetir `run.sh --check` con el backend y las rutas de ese
equipo. Esto verifica la preparación estática; no sustituye una ejecución
limpia ni acredita por sí solo toda la equivalencia científica entre entornos.

# METATEST — release reproducible para CEDIA

Este repositorio contiene el código, la configuración y los manifiestos del
release `metatest-reproducible-20260830`. Los datos, resultados y dos artefactos
sellados de gran tamaño se distribuyen fuera de Git.

## Estado científico y de cierre

- Referencia científica vigente: job **24178**, escenario
  `P1_ALL_qtx00_qpr00_qme00_qcl00` (`q00`), sobre el commit
  `206ac55df4d9a017293b8e7de9f17ca06fbb883a`.
- El job 24178 terminó la validación de q00 mediante una reanudación controlada;
  no constituye por sí solo una prueba desde resultados vacíos.
- El sellado fue validado por el job 24225: imagen, backend MOFA de solo lectura
  y estado Basilisk efímero escribible.
- El cierre científico definitivo del release requiere todavía ejecutar una vez
  `q00` desde un clon limpio y un directorio de resultados vacío. Hasta que esa
  aceptación termine, no debe presentarse el release como validado de extremo a
  extremo desde cero.

## Contenido versionado y privado

Git incluye el código R/shell, plantillas, `renv.lock`, configuración y
manifiestos. No incluye:

- los seis archivos de datos de entrada;
- resultados del pipeline;
- la imagen SquashFS de 5.83 GiB;
- el backend MOFA de 2.02 GiB;
- credenciales ni rutas personales.

El acceso a los datos debe gestionarse por el canal autorizado del proyecto.
Los hashes permiten verificar una copia ya autorizada; no conceden acceso ni
reemplazan las obligaciones de protección de datos.

## Artefactos sellados

| Artefacto | SHA-256 | Bytes |
|---|---|---:|
| `metatest-runtime-metatest-reproducible-20260830.sqsh` | `87b7206e95ca9fd824703cfdff32a0428c75b7b0c2cf2c80abf75d6550c9dc13` | 6260310016 |
| `metatest-mofa-backend-metatest-reproducible-20260830.tar` | `587657369e8687709f65d0e0c6dcc4a81bfaafeaa9ae022d9f1ac7835c678b22` | 2166446080 |

La imagen se inicia con `enroot create` seguido de `enroot start`. Esto evita
depender de `squashfuse`. El login y el host SLURM tampoco necesitan `Rscript`:
R se ejecuta exclusivamente dentro de la imagen sellada.

El payload extraído del backend se monta en modo solo lectura. Basilisk recibe
un pequeño árbol administrativo escribible y efímero, con enlaces hacia el
payload; ese estado se elimina al terminar el job.

## Requisitos de CEDIA

En login: `bash`, `git`, `sbatch`, `sha256sum` y utilidades POSIX habituales.
En el nodo de cómputo: Enroot con `enroot create`, `unsquashfs`, espacio temporal
y acceso a los datos/artefactos. No se restauran paquetes, no se recompila
CmdStan y se bloquean descargas mediante guardas offline.

Los scripts son compatibles con el Git antiguo observado en el nodo login: no
dependen de `git -C`, `git worktree`, porcelain v2 ni opciones recientes para
resolver tags.

La asignación observada en la validación previa fue `cpu-dev`, 8 CPUs, 64 GiB y
8 horas. Es el valor por defecto, no una garantía de duración para una corrida
limpia: q00 podría requerir un límite mayor autorizado por CEDIA o una
reanudación explícita y documentada.

## Estructura recomendada

```text
workspace/
├── metatest/                         # clon limpio y tag del release
├── data/                             # seis inputs autorizados
├── results/                          # salidas, logs y evidencia
└── metatest_release_assets_private/
    └── metatest-reproducible-20260830/
        ├── metatest-runtime-metatest-reproducible-20260830.sqsh
        ├── metatest-mofa-backend-metatest-reproducible-20260830.tar
        ├── release_manifest.env
        ├── inputs.sha256
        ├── SHA256SUMS_RELEASE.txt
        ├── renv.lock.validated
        └── renv.settings.validated.json
```

Las rutas pueden configurarse sin editar el repositorio:

```bash
export METATEST_ASSET_DIR=/ruta/privada/metatest-reproducible-20260830
export METATEST_DATA_DIR=/ruta/autorizada/data
export METATEST_RESULTS_DIR=/ruta/results
```

## Ejecución

Desde un clon limpio situado en el tag `metatest-reproducible-20260830`:

```bash
bash scripts/reproducibility/submit_cedia.sh q00
```

El modo `q00` filtra únicamente el escenario validado. Crea un directorio de
corrida único y exige que `pipeline_results/` esté vacío. Dentro del job se
verifican los hashes completos de imagen, backend e inputs; las identidades de
R, renv, Bioconductor, MOFA2, Basilisk, reticulate, mofapy2, CmdStan,
cmdstanr y ggplot2; las políticas de montaje; el preflight; ambos frameworks;
los heatmaps o su ausencia justificada; la columna `selected_query`; y
`resumen_corridas.tsv`.

La interfaz para todas las filas activas es:

```bash
bash scripts/reproducibility/submit_cedia.sh all
```

Ese modo no forma parte de la aceptación q00 y no debe lanzarse durante este
cierre. Ejecuta las 201 filas activas (67 por cada una de las tres plantillas)
sin cambiar metodología ni parámetros estadísticos.

Para reanudar deliberadamente una corrida existente:

```bash
export METATEST_RESUME=1
export METATEST_RUN_DIR=/ruta/results/reproducibilidad_release/metatest-reproducible-20260830/corrida_previa
bash scripts/reproducibility/submit_cedia.sh q00
```

No se reutiliza un directorio no vacío salvo con esas dos variables. La
aceptación final limpia debe ejecutarse sin ellas.

Los recursos SLURM pueden ajustarse con
`METATEST_SLURM_PARTITION`, `METATEST_SLURM_CPUS`,
`METATEST_SLURM_MEM` y `METATEST_SLURM_TIME`. Los directorios de cada corrida
conservan stdout/stderr, tiempos, estado final, consumo observado y artefactos
de validación.

## Identidades principales

| Componente | Versión/identidad |
|---|---|
| R | 4.4.2 |
| renv | 1.2.0 |
| Bioconductor | 3.20 |
| MOFA2 / Basilisk / reticulate | 1.16.0 / 1.18.0 / 1.45.0 |
| mofapy2 | 0.7.0 |
| CmdStan / cmdstanr | 2.36.0 / 0.8.0 |
| commit remoto de cmdstanr | `12fe0f8fd35226dadf3abc50de726537dc9ea475` |
| ggplot2 | 3.5.2 |
| lock efectivo | `5af95568c692717baa2eb251a699fc706a2de3cf7e74c46441d096435c813b30` |

Las identidades completas y los hashes de los seis inputs están en
`config/cedia_release.env` y `reproducibilidad/release/`.

## Límites metodológicos

Este release no cambia plantillas, parámetros estadísticos ni metodología. En
particular, no modifica `crear_mapa_clinico_maestro.R`,
`crear_mapa_proteico_maestro.R` ni relaja la validación del framework. Los mapas
de referencia se regeneran como salidas ignoradas por Git; por eso el código se
monta escribible, mientras datos y backend permanecen en solo lectura.

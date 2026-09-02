#!/usr/bin/env bash
# Funciones exclusivas del backend nativo. No instalar ni restaurar dependencias.
metatest_native_error() { printf 'ERROR NATIVO: %s\n' "$*" >&2; return 1; }

metatest_native_prepare() {
  local repo="$1" report="$2" cmd candidate count=0 selected='' key value records
  for cmd in Rscript python3 g++ make sha256sum realpath flock; do
    command -v "$cmd" >/dev/null || { metatest_native_error "falta $cmd"; return 1; }
  done
  [[ "$(type -t Rscript)" == file ]] || { metatest_native_error 'Rscript debe ser un ejecutable, no una funcion'; return 1; }
  METATEST_RSCRIPT="$(realpath -e "$(command -v Rscript)")" || return 1
  export METATEST_RSCRIPT
  export METATEST_REPO_ROOT="$repo"
  if [[ -n "${METATEST_R_LIBRARY:-}" ]]; then
    selected="$METATEST_R_LIBRARY"
  else
    # No buscar en bibliotecas globales ni activar renv/activate.R.
    for candidate in "$repo"/renv/library/*/R-4.4 "$repo"/renv/library/*/R-4.4/*; do
      [[ -f "$candidate/renv/DESCRIPTION" && -f "$candidate/jsonlite/DESCRIPTION" ]] || continue
      selected="$candidate"
      count=$((count+1))
    done
    [[ "$count" == 1 ]] || {
      metatest_native_error "se encontraron $count bibliotecas de proyecto R-4.4; restaura renv.lock o indica METATEST_R_LIBRARY explicitamente"
      return 1
    }
  fi
  METATEST_R_LIBRARY="$(realpath -e "$selected")" || return 1
  export METATEST_R_LIBRARY
  [[ -f "$METATEST_R_LIBRARY/renv/DESCRIPTION" ]] || { metatest_native_error 'biblioteca incompleta'; return 1; }
  CMDSTAN="$(realpath -e "${METATEST_CMDSTAN:-${CMDSTAN:-${HOME:?}/.cmdstan/cmdstan-2.36.0}}")" || return 1
  export CMDSTAN
  [[ -d "$CMDSTAN" && -x "$CMDSTAN/bin/stanc" ]] || { metatest_native_error 'falta CmdStan 2.36.0 con bin/stanc'; return 1; }
  for candidate in "$repo" "$report" "$METATEST_R_LIBRARY" "$CMDSTAN"; do
    [[ "$candidate" != *$'\n'* && "$candidate" != *$'\t'* && "$candidate" != *:* ]] || {
      metatest_native_error 'las rutas del runtime no pueden contener TAB, newline ni dos puntos'; return 1;
    }
  done
  # R_LIBS_* tambien se hereda en Rscript --vanilla y en workers nuevos.
  export R_LIBS="$METATEST_R_LIBRARY" R_LIBS_USER="$METATEST_R_LIBRARY" R_LIBS_SITE=NULL
  export R_PROFILE=/dev/null R_PROFILE_USER=/dev/null R_ENVIRON=/dev/null R_ENVIRON_USER=/dev/null
  export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE RENV_CONFIG_AUTO_SNAPSHOT=FALSE
  export METATEST_NATIVE_ENV_ACTIVE=1 CMDSTAN_PATH="$CMDSTAN"
  export PIPELINE_R_CORES=2 BRMS_INTERNAL_CORES=1 BLAS_THREADS=1
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 STAN_NUM_THREADS=1 TBB_NUM_THREADS=1
  export MC_CORES=1 O=0 MAKEFLAGS=-j1 CMDSTAN_OPTIMIZATION_LEVEL=0
  export BASILISK_USE_SYSTEM_DIR=0 BASILISK_NO_DESTROY=1
  export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 RETICULATE_USE_MANAGED_VENV=no
  unset PYTHONPATH PYTHONHOME VIRTUAL_ENV RETICULATE_PYTHON_ENV RETICULATE_PYTHON_FALLBACK
  mkdir -p "$report" || return 1
  records="$report/native_paths.tsv"
  "$METATEST_RSCRIPT" --vanilla "$repo/scripts/portable/native_check.R" resolve > "$records" || return 1
  count=0
  while IFS=$'\t' read -r key value; do
    [[ -n "$value" && "$value" == /* && "$value" != *$'\t'* && "$value" != *:* ]] || {
      metatest_native_error 'salida de resolucion de rutas invalida'; return 1;
    }
    case "$key" in
      BASILISK_EXTERNAL_DIR) export BASILISK_EXTERNAL_DIR="$value";;
      BASILISK_EXTERNAL_CONDA) export BASILISK_EXTERNAL_CONDA="$value";;
      RETICULATE_PYTHON) export RETICULATE_PYTHON="$value";;
      *) metatest_native_error "clave inesperada: $key"; return 1;;
    esac
    count=$((count+1))
  done < "$records"
  [[ "$count" == 3 && -x "$RETICULATE_PYTHON" && -x "$BASILISK_EXTERNAL_CONDA/bin/conda" ]] || {
    metatest_native_error 'falta el entorno MOFA o su Conda existente; no se instala automaticamente'; return 1;
  }
  printf 'R_LIBRARY=%s\nCMDSTAN=%s\nRETICULATE_PYTHON=%s\nHOST_RSCRIPT_REQUIRED=TRUE\n' \
    "$METATEST_R_LIBRARY" "$CMDSTAN" "$RETICULATE_PYTHON"
}

metatest_native_smoke() {
  local repo="$1" report="$2"
  "$METATEST_RSCRIPT" --vanilla "$repo/scripts/portable/native_check.R" check "$report" || return 1
  # Consulta de flags: no compila, no ejecuta make clean ni escribe make/local.
  make --no-print-directory -s -C "$CMDSTAN" print-compiler-flags > "$report/cmdstan_compiler_flags.txt" || return 1
  [[ "$(awk '/O \(Optimization Level\)/ {print $NF; exit}' "$report/cmdstan_compiler_flags.txt")" == 0 ]] || {
    metatest_native_error 'CmdStan no acepta O=0 en este perfil'; return 1;
  }
  g++ --version > "$report/cxx_version.txt" || return 1
  (
    cd "$report" || exit 1
    sha256sum native_identity.json cmdstan_compiler_flags.txt cxx_version.txt
  ) > "$report/native_identity_parts.sha256" || return 1
  METATEST_NATIVE_FINGERPRINT="$(sha256sum "$report/native_identity_parts.sha256" | awk '{print $1}')" || return 1
  export METATEST_NATIVE_FINGERPRINT
  printf 'NATIVE_RUNTIME_FINGERPRINT=%s\nNATIVE_RUNTIME_SMOKE=PASS\nSTAN_COMPILATION_TEST=NOT_RUN\nPIPELINE_EXECUTED=FALSE\n' \
    "$METATEST_NATIVE_FINGERPRINT"
}

metatest_native_guard_resume() {
  local state="$1" resumed="$2" identity="$1/native_runtime.sha256"
  if [[ -e "$identity" ]]; then
    [[ "$(< "$identity")" == "$METATEST_NATIVE_FINGERPRINT" ]] || {
      metatest_native_error 'el entorno nativo cambio; no se permite resume. Usa otro outdir'; return 1;
    }
  elif [[ "$resumed" == 1 ]]; then
    metatest_native_error 'la corrida anterior no tiene huella nativa; usa otro outdir'; return 1
  else
    printf '%s\n' "$METATEST_NATIVE_FINGERPRINT" > "$identity"
  fi
}

#!/usr/bin/env bash
# Verificacion del checkout de trabajo. El runtime y el nucleo cientifico
# conservan la identidad del release; HEAD puede avanzar para trabajo post hoc.
metatest_verify_checkout() {
  local repo="$1" base=3dc20d7faf86d2341988470b20b748cb0fac4255
  local tag=metatest-reproducible-20260830 head ancestor path
  head="$(cd "$repo" && git rev-parse HEAD)" || return 1
  [[ "$(cd "$repo" && git rev-parse "${tag}^{commit}")" == "$base" ]] || {
    printf 'ERROR: identidad del tag historico distinta\n' >&2; return 1;
  }
  ancestor="$(cd "$repo" && git merge-base "$base" "$head")" || return 1
  [[ "$ancestor" == "$base" ]] || {
    printf 'ERROR: HEAD no desciende del release validado\n' >&2; return 1;
  }
  (cd "$repo" && git diff --no-ext-diff --no-textconv --quiet HEAD --) || {
    printf 'ERROR: hay cambios tracked; hacer commit antes de ejecutar\n' >&2; return 1;
  }
  # Las nuevas tareas post hoc pueden ir en scripts/selection u otro directorio.
  # Cambiar estos componentes requiere una nueva validacion, no quitar el guard.
  (cd "$repo" && git diff --no-ext-diff --no-textconv --quiet "$base" HEAD -- \
    scripts/pipeline scripts/utilities templates/current renv.lock renv/settings.json \
    config/cedia_release.env reproducibilidad/release \
    scripts/reproducibility/01_preflight_servidor.sh) || {
    printf 'ERROR: nucleo cientifico/runtime distinto del release; requiere nueva validacion\n' >&2
    return 1
  }
  for path in run.sh run.sbatch PORTABLE_SHA256SUMS.txt scripts/portable/provenance.sh \
    scripts/portable/enroot.sh scripts/portable/local.sh scripts/portable/orchestrator.sh \
    scripts/portable/validate.sh scripts/portable/collect.sh; do
    (cd "$repo" && git ls-files --error-unmatch -- "$path" >/dev/null) || {
      printf 'ERROR: archivo portable sin versionar: %s\n' "$path" >&2; return 1;
    }
  done
  (cd "$repo" && sha256sum -c PORTABLE_SHA256SUMS.txt) || return 1
  export METATEST_CODE_COMMIT="$head"
  export METATEST_RUNTIME_BASE_COMMIT="$base"
  export METATEST_CORE_RELEASE_MATCH=PASS
  export METATEST_PORTABLE_MANIFEST_SHA256="$(sha256sum "$repo/PORTABLE_SHA256SUMS.txt" | awk '{print $1}')"
}

metatest_print_provenance() {
  printf 'CODE_COMMIT=%s\nRUNTIME_BASE_COMMIT=%s\nCORE_MATCHES_RELEASE=%s\nPORTABLE_MANIFEST_SHA256=%s\n' \
    "${METATEST_CODE_COMMIT:?}" "${METATEST_RUNTIME_BASE_COMMIT:?}" \
    "${METATEST_CORE_RELEASE_MATCH:?}" "${METATEST_PORTABLE_MANIFEST_SHA256:?}"
  printf 'SCIENTIFIC_REFERENCE_JOB=24178\nCLEAN_Q00_REFERENCE_JOB=24294\n'
}

#!/usr/bin/env Rscript
# METATEST2 / script 30 EXTENDIDO — CORE28/SHARED/SPECIFIC/ALL + CORE18_47
# y auditoria de pesos MOFA nativos de los modelos ganadores.
# No reentrena MOFA, no recalcula Q5/PLSDA/ElasticNet/SGL ni reejecuta el 28.
# Datos: exportacion COMPLETA del script 28 y ready_for_modeling.rds existentes.
# Todos los calculos predictivos son nuevos y condicionados a la seleccion previa.
# Validacion: CV anidada exploratoria dentro de TRAIN y TEST descriptivo/reutilizado.

options(stringsAsFactors = FALSE)
CLASS_LEVELS <- c("NP", "FD", "SO")
FAMILIES <- c("P1", "P2", "P3")
VIEWS <- c("tx", "pr", "me", "cl")
PANEL_NAMES <- c("CORE", "SHARED", "SPECIFIC", "ALL", "CORE18_47")
# Si ya ejecutaste los 12 modelos originales, ejecutar SOLO los 3 del panel CORE18_47:
# MT2_30_ONLY_CORE18_47=1 Rscript --vanilla <este_script.R>
ONLY_CORE18_47 <- identical(Sys.getenv("MT2_30_ONLY_CORE18_47", "0"), "1")
RUN_PANELS <- if (ONLY_CORE18_47) "CORE18_47" else PANEL_NAMES
N_CV_REPEATS <- as.integer(Sys.getenv("MT2_30_CV_REPEATS", "5"))
N_OUTER_REPEATS <- as.integer(Sys.getenv("MT2_30_OUTER_REPEATS", "3"))
SEED <- 20261101L
stopifnot(N_CV_REPEATS >= 2, N_OUTER_REPEATS >= 1)

# El script vive en <repo>/scripts/biomarker_selection y detecta la raiz de Git
# a partir de su propia ubicacion (independiente del directorio de trabajo).
SCRIPT30_VERSION <- "1.1.0-repo"
script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file_arg)) {
  script_file <- normalizePath(sub("^--file=", "", script_file_arg[1L]), mustWork = TRUE)
  repo_root <- dirname(dirname(dirname(script_file)))
} else {
  script_file <- NA_character_
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
RUN <- path.expand(Sys.getenv("METATEST2_RUN_ROOT", file.path(
  dirname(repo_root), "results_cedia", "metatest2_integral_66_20261027T013422Z")))
SELECTION_RUN <- path.expand(Sys.getenv("METATEST2_SELECTION_RUN", file.path(
  RUN, "comparacion_modelamientos", "upstream_comun", "downstream",
  "seleccion_mejor_escenario_train", "20261029T210729Z_pid2405301")))
assert <- function(cond, ...) if (!isTRUE(cond)) stop(paste(..., collapse=""), call.=FALSE)
rebase <- function(p) {
  p <- as.character(p)[1L]
  assert(!is.na(p) && nzchar(p), "Ruta de entrada vacia")
  prefix <- "/workspace/results/"
  if (startsWith(p,prefix)) return(file.path(RUN, substring(p,nchar(prefix)+1L)))
  path.expand(p)
}
read_csv <- function(p, cols=character()) {
  assert(file.exists(p), "Archivo ausente: ", p)
  d <- utils::read.csv(p, check.names=FALSE, stringsAsFactors=FALSE)
  miss <- setdiff(cols, names(d))
  assert(length(miss)==0L, "Columnas ausentes en ",p,": ",paste(miss,collapse=", "))
  d
}
write_csv <- function(d,p) utils::write.csv(d,p,row.names=FALSE, na="", fileEncoding="UTF-8")
latest_dir <- function(root, predicate) {
  assert(dir.exists(root), "Directorio ausente: ",root)
  dirs <- list.dirs(root, recursive=FALSE, full.names=TRUE)
  dirs <- dirs[vapply(dirs, predicate, logical(1))]
  assert(length(dirs)>0L,"Ninguna corrida valida dentro de ",root)
  dirs[order(file.info(dirs)$mtime, basename(dirs), decreasing=TRUE)][1L]
}

export_override <- Sys.getenv("METATEST2_SCRIPT28_EXPORT_DIR", "")
if (nzchar(export_override)) {
  export_root <- path.expand(export_override)
} else {
  export_root <- latest_dir(file.path(SELECTION_RUN,"FINAL_STORY_FIGURES_TABLES"),
    function(d) all(file.exists(file.path(d,"00_GLOBAL","tables",c(
      "FINAL_FEATURES_ALL.csv", "FINAL_FEATURE_SHARING.csv",
      "MOFA_NATIVE_IMPORTANCE.csv", "EVIDENCE_A_ALL_MODEL_VOTES.csv")))))
}
SRC <- file.path(export_root,"00_GLOBAL","tables")
assert(dir.exists(SRC),"Falta exportacion 28: ",SRC)
final <- read_csv(file.path(SRC,"FINAL_FEATURES_ALL.csv"),
  c("family","view_code","feature_original","feature_id","predictor_key"))
sharing <- read_csv(file.path(SRC,"FINAL_FEATURE_SHARING.csv"),
  c("feature_id","n_questions"))
mofa <- read_csv(file.path(SRC,"MOFA_NATIVE_IMPORTANCE.csv"),
  c("family","view_code","feature_original","MOFA_loading_norm"))
evidence <- read_csv(file.path(SRC,"EVIDENCE_A_ALL_MODEL_VOTES.csv"),
  c("family","view_code","feature_original"))
assert(!anyDuplicated(final[c("family","view_code","feature_original")]),
  "Features finales duplicadas por familia, vista e ID")
assert(!anyDuplicated(mofa[c("family","view_code","feature_original")]),
  "MOFA_NATIVE_IMPORTANCE contiene IDs duplicados")
assert(!anyDuplicated(sharing$feature_id),"SHARING contiene IDs duplicados")
assert(all(final$family %in% FAMILIES) && all(final$view_code %in% VIEWS),
  "Familia/vista no reconocida")
final$n_questions <- sharing$n_questions[match(final$feature_id,sharing$feature_id)]
assert(all(final$n_questions %in% 1:3),"Inconsistencia al enlazar tabla sharing")
final$panel_class <- c("SPECIFIC","SHARED","CORE")[final$n_questions]
assert(all(table(final$family) > 0),"Una familia sin features")

# Los nombres y la norma MOFA proceden del 28, no se recalculan pesos.
OUT <- file.path(SELECTION_RUN,"SCRIPT30_LOGISTICA_MOFA",
  paste0(format(Sys.time(),"%Y%m%dT%H%M%S"),"_pid",Sys.getpid()))
dir.create(file.path(OUT,"tables"), recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(OUT,"figures"),recursive=TRUE,showWarnings=FALSE)
DIRT <- file.path(OUT,"tables"); DIRF <- file.path(OUT,"figures")
dir.create(file.path(OUT,"models"),recursive=TRUE,showWarnings=FALSE)
git_commit <- tryCatch({
  if (!dir.exists(file.path(repo_root, ".git"))) NA_character_ else {
    z <- suppressWarnings(system2("git", c("-C", shQuote(repo_root), "rev-parse", "HEAD"),
                                  stdout = TRUE, stderr = FALSE))
    if (length(z) && is.null(attr(z, "status"))) z[[1L]] else NA_character_
  }
}, error = function(e) NA_character_)
writeLines(c(
  paste("SCRIPT30_VERSION:", SCRIPT30_VERSION),
  paste("SCRIPT30_FILE:", script_file),
  paste("SCRIPT30_MD5:", if (is.na(script_file)) NA_character_ else
        unname(tools::md5sum(script_file))),
  paste("GIT_HEAD_AT_START:", git_commit),
  paste("Script 28:",normalizePath(SRC)),
  paste("RUN:",normalizePath(RUN)),
  "MOFA: pesos previamente exportados; NO reentrenamiento",
  "CORE=3/3 del script 28 (33 comunes), SHARED=exactamente 2/3, SPECIFIC=1/3, ALL=union del 28",
  "CORE18_47=47 biomarcadores CORE del script 18, SIN filtrar por votos del script 28",
  paste("Paneles ejecutados:",paste(RUN_PANELS,collapse=", ")),
  "Ridge multinomial (alpha=0), regularizacion CV TRAIN; test reutilizado exploratorio",
  "La seleccion biologica previa al CV condiciona la CV; no validacion anidada de todo el pipeline"
),file.path(OUT,"README.txt"))

# UMBRAL UNICO EN METODO: Otsu sin top-N ni cuantiles subjetivos.
# Se aplica por familia/vista a w^2 porque los pesos de distintas vistas no
# deben compararse en bruto. Un corte Otsu es descriptivo, NO p-value/FDR.
otsu_cut <- function(values) {
  x <- sort(unique(as.numeric(values[is.finite(values)])))
  if (length(x)<2L) return(c(threshold=NA_real_, separation=NA_real_))
  v <- as.numeric(values[is.finite(values)])
  total_var <- mean((v-mean(v))^2)
  candidates <- (x[-length(x)]+x[-1L])/2
  between <- vapply(candidates, function(t) {
    a <- v[v<=t]; b <- v[v>t]
    if (!length(a)||!length(b)) return(NA_real_)
    (length(a)*length(b)/length(v)^2)*(mean(a)-mean(b))^2
  },numeric(1))
  i <- which.max(between)
  c(threshold=candidates[i], separation=if (total_var>0) between[i]/total_var else NA_real_)
}
mo <- mofa
mo$MOFA_loading_norm <- suppressWarnings(as.numeric(mo$MOFA_loading_norm))
assert(all(is.finite(mo$MOFA_loading_norm) & mo$MOFA_loading_norm>=0),
  "Pesos MOFA no finitos o negativos")
mo$key <- paste(mo$family,mo$view_code,mo$feature_original,sep="\r")
final$key <- paste(final$family,final$view_code,final$feature_original,sep="\r")
mo$in_final <- mo$key %in% final$key
mo$panel_class <- final$panel_class[match(mo$key,final$key)]
mo$panel_class[is.na(mo$panel_class)] <- "NOT_FINAL"
mo$above_mofa_threshold <- NA
thresholds <- list(); overlap <- list(); counter <- 1L
for (fam in FAMILIES) for (v in VIEWS) {
  i <- which(mo$family==fam & mo$view_code==v)
  z <- mo[i,,drop=FALSE]
  if (!nrow(z)) next
  cut <- otsu_cut(z$MOFA_loading_norm^2)
  thr <- sqrt(cut["threshold"])
  mo$above_mofa_threshold[i] <- if (is.finite(thr))
    z$MOFA_loading_norm > thr else NA
  thresholds[[length(thresholds)+1L]] <- data.frame(
    family=fam,view_code=v,n_features=nrow(z),otsu_threshold_norm=unname(thr),
    otsu_separation_ratio=unname(cut["separation"]),
    n_above=sum(mo$above_mofa_threshold[i],na.rm=TRUE))
  sel <- mo$in_final[i]; hi <- mo$above_mofa_threshold[i]
  overlap[[length(overlap)+1L]] <- data.frame(
    family=fam,view_code=v,n_model=nrow(z),n_final=sum(sel),
    n_above=sum(hi,na.rm=TRUE),n_final_and_above=sum(sel & hi,na.rm=TRUE),
    n_final_below=sum(sel & !hi,na.rm=TRUE),
    n_above_not_final=sum(!sel & hi,na.rm=TRUE))
  z <- mo[i,,drop=FALSE]
  z <- z[order(-z$MOFA_loading_norm,z$feature_original),,drop=FALSE]
  png(file.path(DIRF,paste0("MOFA_",fam,"_",toupper(v),"_ALL.png")),
      width=2200,height=1250,res=165)
  par(mar=c(if(nrow(z)<=38) 13 else 5, 5, 5, 2))
  cc <- ifelse(z$in_final,"#2370B7","#BBBBBB")
  plot(seq_len(nrow(z)),z$MOFA_loading_norm,col=cc,pch=19,cex=0.85,
       xaxt="n",xlab="Features ordenadas por norma del peso MOFA",
       ylab="Norma L2 de pesos (2 factores)",
       main=paste("MOFA nativo:",fam,toupper(v)))
  abline(h=unname(thr),col="#BB3A2B",lwd=2,lty=2)
  if(nrow(z)<=38) axis(1,at=seq_len(nrow(z)),labels=z$feature_original,las=2,cex.axis=0.65)
  legend("topright",legend=c("Panel 28","No final","Umbral Otsu"),
    col=c("#2370B7","#BBBBBB","#BB3A2B"),pch=c(19,19,NA),
    lty=c(NA,NA,2),bty="n",cex=0.85)
  dev.off()
  selected <- z[z$in_final,,drop=FALSE]
  if(nrow(selected)) {
    selected <- selected[order(selected$MOFA_loading_norm),,drop=FALSE]
    png(file.path(DIRF,paste0("MOFA_",fam,"_",toupper(v),"_FINAL_ONLY.png")),
      width=2200,height=max(850,round(140*nrow(selected))),res=165)
    par(mar=c(5, min(15,max(7,max(nchar(selected$feature_original))/3)),5,2))
    barplot(selected$MOFA_loading_norm,names.arg=selected$feature_original,
      horiz=TRUE,las=1,col=ifelse(selected$above_mofa_threshold,"#2370B7","#B9B9B9"),
      border=NA,cex.names=0.8, xlab="Norma L2 de pesos MOFA",
      main=paste(fam,toupper(v),"— features finales del 28"))
    abline(v=unname(thr),col="#BB3A2B",lwd=2,lty=2)
    dev.off()
  }
}
write_csv(mo[,c("family","view_code","feature_original","MOFA_loading_norm",
               "in_final","panel_class","above_mofa_threshold")],
          file.path(DIRT,"01_MOFA_FEATURE_VS_PANEL28.csv"))
write_csv(do.call(rbind,thresholds),file.path(DIRT,"02_MOFA_THRESHOLDS_OTSU.csv"))
write_csv(do.call(rbind,overlap),file.path(DIRT,"03_MOFA_OVERLAP_BY_VIEW.csv"))

# Rastreo visible de las dos variables clinicas de interes, sin inferir causas.
trace <- do.call(rbind,lapply(FAMILIES,function(fam) do.call(rbind,lapply(c("FLI","PCR"),function(feature) {
  key <- paste(fam,"cl",feature,sep="\r")
  iz <- which(mo$key==key)
  data.frame(family=fam,feature=feature,in_mofa=length(iz)>0,
    mofa_norm=if(length(iz)) mo$MOFA_loading_norm[iz[1L]] else NA_real_,
    above_otsu=if(length(iz)) mo$above_mofa_threshold[iz[1L]] else NA,
    in_evidence_A=any(evidence$family==fam & evidence$view_code=="cl" &
      evidence$feature_original==feature),
    in_final=key %in% final$key)
}))))
write_csv(trace,file.path(DIRT,"04_FLI_PCR_TRACE.csv"))
write_csv(as.data.frame(table(final$family,final$panel_class)),
  file.path(DIRT,"05_PANEL_SIZE_AUDIT.csv"))

# PESOS FIRMADOS POR FACTOR (opcional): leer el RDS MOFA preexistente.
# No hay reentrenamiento; los pesos firmados pueden cambiar de signo por la
# indeterminacion de orientacion de factores. La regla Otsu SOLO se aplica
# a las normas L2 del CSV previo, no se escogen umbrales extra por factor.
if (requireNamespace("MOFA2",quietly=TRUE)) {
  signed_list <- list()
  for (fam in FAMILIES) {
    root <- file.path(SELECTION_RUN,fam,"BIDIRECTIONAL_FRAMEWORK")
    audits <- if (dir.exists(root)) list.files(root,pattern="^00_INPUT_AUDIT\\.csv$",
        recursive=TRUE,full.names=TRUE) else character()
    audits <- audits[grepl("/CROSSOMIC_BRIDGE/.*/tables/00_INPUT_AUDIT\\.csv$",audits)]
    if (!length(audits)) {
      warning("No se encontro CROSSOMIC_BRIDGE audit para pesos firmados: ",fam)
      next
    }
    aud <- read_csv(audits[order(file.info(audits)$mtime,decreasing=TRUE)[1L]],
      "model_path")
    mp <- rebase(aud$model_path[1L])
    if (!file.exists(mp)) {
      warning("No existe RDS de MOFA para pesos firmados: ",mp)
      next
    }
    w <- MOFA2::get_weights(readRDS(mp),views="all",factors="all",as.data.frame=TRUE)
    assert(all(c("view","feature","factor","value") %in% names(w)),
      "MOFA2::get_weights devolvio estructura inesperada")
    map_view <- function(x) {
      y <- tolower(iconv(as.character(x),to="ASCII//TRANSLIT"))
      ifelse(grepl("transcript",y),"tx",ifelse(grepl("prote",y),"pr",
        ifelse(grepl("metabol",y),"me",ifelse(grepl("clinic",y),"cl",NA_character_))))
    }
    w$view_code <- map_view(w$view)
    assert(all(!is.na(w$view_code)),"MOFA2 contiene nombres de vistas no reconocidos")
    w$family <- fam
    w$feature_original <- as.character(w$feature)
    w$value <- as.numeric(w$value)
    w$key <- paste(w$family,w$view_code,w$feature_original,sep="\r")
    w$in_final <- w$key %in% final$key
    w$above_mofa_threshold <- mo$above_mofa_threshold[match(w$key,mo$key)]
    signed_list[[fam]] <- w[,c("family","view_code","feature_original","factor",
                               "value","in_final","above_mofa_threshold")]
    for (v in VIEWS) for (fac in unique(w$factor)) {
      zz <- w[w$view_code==v & w$factor==fac,,drop=FALSE]
      if (!nrow(zz)) next
      zz <- zz[order(-abs(zz$value),zz$feature_original),,drop=FALSE]
      # El tope 30 SOLO controla etiquetas, no filtra datos ni seleccion biologica.
      zz <- head(zz,30L)
      zz <- zz[order(zz$value),,drop=FALSE]
      png(file.path(DIRF,paste0("MOFA_SIGNED_",fam,"_",toupper(v),"_",fac,".png")),
        width=2000,height=1900,res=160)
      par(mar=c(5,11,5,2))
      barplot(zz$value,horiz=TRUE,names.arg=zz$feature_original,las=1,
        col=ifelse(zz$in_final,"#2370B7","#B0B0B0"),border=NA,
        cex.names=0.65,xlab="Peso MOFA firmado (signo arbitrario)",
        main=paste(fam,toupper(v),fac,"— top 30 por |W| (solo vista)"))
      abline(v=0,col="black")
      dev.off()
    }
  }
  if (length(signed_list)) {
    signed_all <- do.call(rbind,signed_list)
    write_csv(signed_all,file.path(DIRT,"06_MOFA_SIGNED_WEIGHTS_BY_FACTOR.csv"))
    # Comprobar que la norma del factor exportado coincide con el paso 28.
    norm <- aggregate(I(value^2)~family+view_code+feature_original,signed_all,sum)
    names(norm)[4L] <- "sum_squared"
    norm$norm_from_factors <- sqrt(norm$sum_squared)
    orig <- merge(norm, mo[,c("family","view_code","feature_original",
        "MOFA_loading_norm")],by=c("family","view_code","feature_original"))
    orig$abs_delta <- abs(orig$norm_from_factors-orig$MOFA_loading_norm)
    write_csv(orig,file.path(DIRT,"07_MOFA_NORM_CONSISTENCY.csv"))
    assert(nrow(orig)>0L && max(orig$abs_delta)<1e-7,
      "Pesos firmados y normas nativas del 28 no coinciden: comprueba modelos y rutas")
  }
} else {
  warning("MOFA2 no instalado: se mantienen graficos L2 y logistica; sin plots firmados")
}

# MODELAMIENTO: reutilizar ready_for_modeling del pipeline, sin calcular matrices
# o imputaciones nuevas sobre datos biologicos originales. El script 10 contiene
# el mismo recorrido de auditorias y orientacion matriz features x muestras.
extract_meta <- function(g) {
  if (is.data.frame(g)) {
    assert(all(c("id","grupo") %in% names(g)),"Metadatos sin id/grupo")
    return(data.frame(id=as.character(g$id),grupo=as.character(g$grupo)))
  }
  assert(is.atomic(g) && !is.null(names(g)),"Estructura de grupos no reconocida")
  data.frame(id=names(g),grupo=as.character(g))
}
load_family <- function(fam,bio_run) {
  override <- Sys.getenv(paste0("METATEST2_RFM_",fam),"")
  if(nzchar(override)) {
    path <- path.expand(override)
  } else {
    audit <- read_csv(file.path(bio_run,fam,"00_PROGRAM_DEFINITION_AUDIT.csv"),
      "observed_train_run")
    obs <- rebase(audit$observed_train_run[1L])
    audit2 <- read_csv(file.path(obs,"tables","00_INPUT_AUDIT.csv"),"ready_for_modeling")
    path <- rebase(audit2$ready_for_modeling[1L])
  }
  assert(file.exists(path),"No existe ready_for_modeling (",fam,"): ",path,
         "\nPuedes establecer METATEST2_RFM_",fam," con su ruta real")
  rfm <- readRDS(path)
  assert(all(c(VIEWS,"grupo_train","grupo_test") %in% names(rfm)),
    "ready_for_modeling sin vistas/grupos esperados: ",path)
  meta <- list(train=extract_meta(rfm$grupo_train),test=extract_meta(rfm$grupo_test))
  for (s in c("train","test")) {
    assert(!anyDuplicated(meta[[s]]$id),"IDs duplicados ",fam," ",s)
    assert(all(meta[[s]]$grupo %in% CLASS_LEVELS),"Grupos inesperados ",fam," ",s)
    for (v in VIEWS) {
      a <- as.matrix(rfm[[v]][[s]])
      assert(!is.null(rownames(a)) && !is.null(colnames(a)),
        "Nombres de features/muestras ausentes: ",fam," ",v," ",s)
      assert(all(meta[[s]]$id %in% colnames(a)),
        "IDs desalineados en ",fam," ",v," ",s)
      assert(!anyDuplicated(rownames(a)),"Features duplicadas en ",fam," ",v," ",s)
      rfm[[v]][[s]] <- a[,meta[[s]]$id,drop=FALSE]
    }
  }
  list(rfm=rfm,meta=meta,rds_path=path)
}

# El paso descriptivo MOFA puede completarse sin glmnet; la regresion si lo requiere.
assert(requireNamespace("glmnet",quietly=TRUE),
  "glmnet no instalado; ejecuta install.packages('glmnet') en R local")
bio_run <- latest_dir(file.path(SELECTION_RUN,"BIOLOGICAL_PROGRAMS"),
  function(d) file.exists(file.path(d,"P1","00_PROGRAM_DEFINITION_AUDIT.csv")))

# CONSENSO ALTERNATIVO DEL SCRIPT 18: 47 features comunes en las tres familias.
# La fuente es el FINAL_CONSENSUS_PANEL de ESTA MISMA CORRIDA, no la del 28.
# No se vuelve a seleccionar, ni se cruza con el filtro 2/3 votos del 28.
consensus18_run <- latest_dir(file.path(bio_run,"FINAL_CONSENSUS_PANEL"),
  function(d) file.exists(file.path(d,"tables","01_FINAL_PANEL_CORE_3OF3.csv")))
consensus18_source <- file.path(consensus18_run,"tables","01_FINAL_PANEL_CORE_3OF3.csv")
core18 <- read_csv(consensus18_source,
  c("view_code","feature_original","n_families","consensus_tier"))
core18$view_code <- tolower(trimws(as.character(core18$view_code)))
core18$feature_original <- as.character(core18$feature_original)
assert(nrow(core18)==47L,"El CORE del script 18 no tiene exactamente 47 features: ",
       nrow(core18),". No se sustituye por otro consenso")
assert(!anyDuplicated(core18[c("view_code","feature_original")]),
       "CORE18_47 con features duplicadas")
assert(all(core18$n_families==3L & core18$consensus_tier=="CORE"),
       "El panel del 18 no es el consenso CORE exacto 3/3")
assert(all(core18$view_code %in% VIEWS),"Vistas inesperadas en el CORE18_47")
# Validar contra el universo Evidence-A ya cargado, NO contra las features finales del 28.
assert(all(c("family","predictor_key") %in% names(evidence)),
       "Evidence-A no incluye family/predictor_key para identificar CORE18")
core18_key <- paste(core18$view_code,core18$feature_original,sep="\r")
core18_by_family <- list()
core18_audit <- list()
for (fam in FAMILIES) {
  e <- evidence[evidence$family==fam,,drop=FALSE]
  ekey <- paste(as.character(e$view_code),as.character(e$feature_original),sep="\r")
  assert(!anyDuplicated(ekey),"IDs duplicados Evidence-A: ",fam)
  ix <- match(core18_key, ekey)
  assert(!anyNA(ix),"No se pudieron recuperar los 47 predictores en Evidence-A para ",fam,
         "; faltan: ",paste(core18$view_code[is.na(ix)],core18$feature_original[is.na(ix)],
                          sep="::",collapse=", "))
  ff <- core18[,c("view_code","feature_original"),drop=FALSE]
  ff$family <- fam
  ff$predictor_key <- as.character(e$predictor_key[ix])
  assert(!anyNA(ff$predictor_key) & all(nzchar(ff$predictor_key)),
         "predictor_key invalido CORE18_47: ",fam)
  assert(!anyDuplicated(ff$predictor_key),"predictor_key duplicado CORE18_47: ",fam)
  ff$in_final28 <- paste(fam,ff$view_code,ff$feature_original,sep="\r") %in% final$key
  core18_by_family[[fam]] <- ff
  core18_audit[[fam]] <- data.frame(
    family=fam,panel="CORE18_47",view_code=ff$view_code,
    feature_original=ff$feature_original,predictor_key=ff$predictor_key,
    in_final28=ff$in_final28,source18=consensus18_source,
    stringsAsFactors=FALSE)
}
write_csv(do.call(rbind,core18_audit),
          file.path(DIRT,"20_CORE18_47_SOURCE_AND_PANEL28_OVERLAP.csv"))
writeLines(c(paste("Consenso del script 18:",consensus18_source),
             "47 features iguales en P1/P2/P3; 14 TX, 4 PR, 28 ME, 1 CL (PCR).",
             "No intersectar con el CORE del script 28: 33 features, regla >=2/3 votos."),
           file.path(OUT,"CONSENSO_47_ORIGEN.txt"))

# Respetar especificacion del script 10: pesos inversos de frecuencia de clase.
class_weights <- function(y) {
  y <- factor(y,levels=CLASS_LEVELS); tab <- table(y)
  assert(all(tab>0),"La particion TRAIN carece de una clase")
  as.numeric((length(y)/(length(CLASS_LEVELS)*as.numeric(tab)))[match(as.character(y),names(tab))])
}
stratified_folds <- function(y,k,seed) {
  set.seed(seed); out <- integer(length(y))
  for (cls in CLASS_LEVELS) {
    ids <- sample(which(as.character(y)==cls))
    out[ids] <- rep(seq_len(k),length.out=length(ids))
  }
  out
}
matrices_for_panel <- function(dat,panel_features,fam,panel) {
  metas <- dat$meta
  rows <- list(); tr <- list(); te <- list()
  for(j in seq_len(nrow(panel_features))) {
    v <- panel_features$view_code[j]; id <- panel_features$feature_original[j]
    key <- panel_features$predictor_key[j]
    a <- dat$rfm[[v]]$train; b <- dat$rfm[[v]]$test
    reason <- "ELIGIBLE"
    if (!id %in% rownames(a)) reason <- "ABSENT_TRAIN"
    else {
      x <- suppressWarnings(as.numeric(a[id,,drop=TRUE]))
      if(any(!is.finite(x))) reason <- "NONFINITE_TRAIN"
      else if(!is.finite(sd(x)) || sd(x)<=0) reason <- "ZERO_VARIANCE_TRAIN"
    }
    if(reason=="ELIGIBLE") {
      xx <- if(id %in% rownames(b)) suppressWarnings(as.numeric(b[id,,drop=TRUE])) else
          rep(NA_real_,nrow(metas$test))
      xx[!is.finite(xx)] <- stats::median(x)
      tr[[key]] <- x; te[[key]] <- xx
    }
    rows[[j]] <- data.frame(family=fam,panel=panel,predictor_key=key,
      view_code=v,feature_original=id,reason=reason)
  }
  eligibility <- do.call(rbind,rows)
  if(!length(tr)) return(list(eligibility=eligibility,skip="Sin predictores elegibles"))
  xtr <- do.call(cbind,tr); xte <- do.call(cbind,te)
  storage.mode(xtr) <- "double"; storage.mode(xte) <- "double"
  rownames(xtr) <- metas$train$id; rownames(xte) <- metas$test$id
  colnames(xtr) <- names(tr); colnames(xte) <- names(te)
  list(eligibility=eligibility,xtr=xtr,xte=xte,
       ytr=factor(metas$train$grupo,levels=CLASS_LEVELS),
       yte=factor(metas$test$grupo,levels=CLASS_LEVELS))
}

# Full TRAIN: lambda.1se repetida en folds estratificados; una sola alpha=0.
fit_ridge <- function(x,y,seed,nrep=N_CV_REPEATS) {
  k <- min(3L, min(table(y)))
  assert(k>=2L,"Faltan muestras por clase para CV")
  cv <- lapply(seq_len(nrep),function(i) {
    foldid <- stratified_folds(y,k,seed+i)
    glmnet::cv.glmnet(x=x,y=y,family="multinomial",
      type.multinomial="grouped",alpha=0,standardize=TRUE,
      weights=class_weights(y), foldid=foldid, nfolds=k,
      type.measure="deviance")
  })
  lambda <- stats::median(vapply(cv,function(z)z$lambda.1se,numeric(1)))
  fit <- glmnet::glmnet(x=x,y=y,family="multinomial",
      type.multinomial="grouped",alpha=0,standardize=TRUE,
      weights=class_weights(y))
  list(fit=fit,lambda=lambda,lambda_1se=vapply(cv,function(z)z$lambda.1se,numeric(1)),
       cv_deviance=vapply(cv,function(z)z$cvm[which.min(abs(log(z$lambda)-log(z$lambda.1se)))],numeric(1)))
}
prob_matrix <- function(fit,x,lambda) {
  raw <- predict(fit,newx=x,s=lambda,type="response")
  assert(length(dim(raw))==3L,"glmnet no devolvio probabilidades 3D")
  mat <- matrix(raw[,,1L],nrow=nrow(x),ncol=dim(raw)[2L],
      dimnames=list(rownames(x),dimnames(raw)[[2L]]))
  assert(all(CLASS_LEVELS %in% colnames(mat)),"Probabilidades sin clases NP FD SO")
  mat[,CLASS_LEVELS,drop=FALSE]
}
assess <- function(probs,truth,fam,panel,split) {
  y <- factor(truth,levels=CLASS_LEVELS)
  predicted <- factor(CLASS_LEVELS[max.col(probs,ties.method="first")],levels=CLASS_LEVELS)
  conf <- table(truth=y,predicted=predicted)
  rec <- diag(conf)/rowSums(conf); rec[!is.finite(rec)] <- NA_real_
  prec <- diag(conf)/colSums(conf); prec[!is.finite(prec)] <- NA_real_
  f1 <- ifelse((prec+rec)>0,2*prec*rec/(prec+rec),0)
  f1[!is.finite(f1)] <- 0
  obs_prob <- probs[cbind(seq_along(y),match(as.character(y),colnames(probs)))]
  metr <- data.frame(family=fam,panel=panel,split=split,n=length(y),
    accuracy=mean(as.character(predicted)==as.character(y)),
    balanced_accuracy=mean(rec,na.rm=TRUE),macro_f1=mean(f1),
    log_loss=-mean(log(pmax(obs_prob,1e-15))),
    brier=mean(rowSums((probs - outer(as.character(y),CLASS_LEVELS,"=="))^2)))
  classes <- data.frame(family=fam,panel=panel,split=split,class=CLASS_LEVELS,
    n=as.integer(table(y)),recall=as.numeric(rec),precision=as.numeric(prec),f1=as.numeric(f1))
  pred <- data.frame(family=fam,panel=panel,split=split,sample_id=rownames(probs),
    true=as.character(y),predicted=as.character(predicted),
    prob_NP=probs[,"NP"],prob_FD=probs[,"FD"],prob_SO=probs[,"SO"])
  list(metrics=metr,classes=classes,pred=pred,conf=conf)
}

metrics <- list(); class_metrics <- list(); predictions <- list(); confusions <- list()
coefficients <- list(); eligibility <- list(); configs <- list(); status <- list()
for(fam in FAMILIES) {
  message("Leyendo matrices ya preprocesadas: ",fam)
  dat <- load_family(fam,bio_run)
  sz <- c(TRAIN=nrow(dat$meta$train),TEST=nrow(dat$meta$test))
  for(panel in RUN_PANELS) {
    if (panel=="CORE18_47") {
      ff <- core18_by_family[[fam]]
    } else {
      ff <- final[final$family==fam,,drop=FALSE]
      if(panel=="CORE") ff <- ff[ff$n_questions==3L,,drop=FALSE]
      if(panel=="SHARED") ff <- ff[ff$n_questions==2L,,drop=FALSE]
      if(panel=="SPECIFIC") ff <- ff[ff$n_questions==1L,,drop=FALSE]
    }
    key <- paste(fam,panel,sep="_")
    message("Modelo ",key," (",nrow(ff)," features de panel)")
    if(!nrow(ff)) {
      status[[key]] <- data.frame(family=fam,panel=panel,status="SKIP_EMPTY_PANEL",detail="Sin features")
      next
    }
    dd <- matrices_for_panel(dat,ff,fam,panel)
    eligibility[[key]] <- dd$eligibility
    if(!is.null(dd$skip)) {
      status[[key]] <- data.frame(family=fam,panel=panel,status="SKIP_NO_FEATURES",detail=dd$skip)
      next
    }
    # Evaluacion CV anidada condicional, siempre con folds TRAIN:
    # Outer 3-fold, inner 2-fold cuando la clase FD queda con 2 observaciones.
    oof_rows <- list(); rr <- 1L; err <- NULL
    for(rep in seq_len(N_OUTER_REPEATS)) {
      foldid <- stratified_folds(dd$ytr,3L,SEED+10000L*match(fam,FAMILIES)+
                                      1000L*match(panel,PANEL_NAMES)+rep)
      pred_rep <- matrix(NA_real_,nrow(dd$xtr),length(CLASS_LEVELS),
                         dimnames=list(rownames(dd$xtr),CLASS_LEVELS))
      for(k in 1:3) {
        itr <- which(foldid!=k); ite <- which(foldid==k)
        tryCatch({
          inner <- fit_ridge(dd$xtr[itr,,drop=FALSE],dd$ytr[itr],
            SEED+100000L*rep+1000L*match(fam,FAMILIES)+100L*match(panel,PANEL_NAMES)+k,
            nrep=max(2L,ceiling(N_CV_REPEATS/2)))
          pred_rep[ite,] <- prob_matrix(inner$fit,dd$xtr[ite,,drop=FALSE],inner$lambda)
        },error=function(e){err <<- conditionMessage(e)})
        if(!is.null(err)) break
      }
      if(!is.null(err)) break
      stopifnot(all(is.finite(pred_rep)))
      cv_res <- assess(pred_rep,dd$ytr,fam,panel,paste0("TRAIN_NESTED_CV_REPEAT_",rep))
      metrics[[paste0(key,"_CV_",rep)]] <- cv_res$metrics
      class_metrics[[paste0(key,"_CV_",rep)]] <- cv_res$classes
      predictions[[paste0(key,"_CV_",rep)]] <- cv_res$pred
      confusions[[paste0(key,"_CV_",rep)]] <- data.frame(family=fam,panel=panel,
        split=paste0("TRAIN_NESTED_CV_REPEAT_",rep),as.data.frame(cv_res$conf))
    }
    if(!is.null(err)) {
      status[[key]] <- data.frame(family=fam,panel=panel,status="FAIL_NESTED_CV",detail=err)
      next
    }
    fitted <- tryCatch(fit_ridge(dd$xtr,dd$ytr,SEED+match(fam,FAMILIES)*100+match(panel,PANEL_NAMES)),
                       error=function(e)e)
    if(inherits(fitted,"error")) {
      status[[key]] <- data.frame(family=fam,panel=panel,status="FAIL_FINAL_FIT",detail=conditionMessage(fitted))
      next
    }
    test_p <- prob_matrix(fitted$fit,dd$xte,fitted$lambda)
    test_res <- assess(test_p,dd$yte,fam,panel,"TEST_REUSED_EXPLORATORY")
    metrics[[paste0(key,"_TEST")]] <- test_res$metrics
    class_metrics[[paste0(key,"_TEST")]] <- test_res$classes
    predictions[[paste0(key,"_TEST")]] <- test_res$pred
    confusions[[paste0(key,"_TEST")]] <- data.frame(family=fam,panel=panel,
      split="TEST_REUSED_EXPLORATORY",as.data.frame(test_res$conf))
    saveRDS(list(model=fitted$fit,lambda=fitted$lambda,
      predictor_keys=colnames(dd$xtr),family=fam,panel=panel,
      n_train=nrow(dd$xtr),classes=CLASS_LEVELS),
      file.path(OUT,"models",paste0(key,"_ridge_multinomial.rds")))
    cf <- stats::coef(fitted$fit,s=fitted$lambda)
    coefficients[[key]] <- do.call(rbind,lapply(names(cf),function(cls) {
      cc <- as.matrix(cf[[cls]])
      data.frame(family=fam,panel=panel,class=cls,predictor_key=rownames(cc),
                 coefficient=as.numeric(cc[,1L]))
    }))
    configs[[key]] <- data.frame(family=fam,panel=panel,n_train=nrow(dd$xtr),
      n_test=nrow(dd$xte),n_selected=nrow(ff),n_eligible=ncol(dd$xtr),
      lambda=fitted$lambda,lambda_1se_median=median(fitted$lambda_1se),
      cv_deviance_median=median(fitted$cv_deviance),
      alpha=0,n_cv_repeats=N_CV_REPEATS,n_outer_repeats=N_OUTER_REPEATS,
      ready_for_modeling=dat$rds_path)
    status[[key]] <- data.frame(family=fam,panel=panel,status="PASS",detail="Final TRAIN y TEST exploratorio")
  }
}

bind_or_empty <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
write_csv(bind_or_empty(status),file.path(DIRT,"10_MODEL_STATUS.csv"))
write_csv(bind_or_empty(eligibility),file.path(DIRT,"11_MODEL_FEATURE_ELIGIBILITY.csv"))
write_csv(bind_or_empty(configs),file.path(DIRT,"12_MODEL_HYPERPARAMETERS.csv"))
write_csv(bind_or_empty(metrics),file.path(DIRT,"13_MODEL_METRICS.csv"))
write_csv(bind_or_empty(class_metrics),file.path(DIRT,"14_PER_CLASS_METRICS.csv"))
write_csv(bind_or_empty(predictions),file.path(DIRT,"15_PREDICTIONS_PER_SAMPLE.csv"))
write_csv(bind_or_empty(coefficients),file.path(DIRT,"16_LOGISTIC_COEFFICIENTS.csv"))
write_csv(bind_or_empty(confusions),file.path(DIRT,"19_CONFUSION_LONG.csv"))

if(length(metrics)) {
  mm <- bind_or_empty(metrics)
  for(split_regex in c("^TRAIN_NESTED_CV","^TEST_REUSED")) {
    ss <- mm[grepl(split_regex,mm$split),,drop=FALSE]
    if(!nrow(ss)) next
    ag <- stats::aggregate(cbind(accuracy,balanced_accuracy,macro_f1) ~ family + panel,
      ss,FUN=mean)
    write_csv(ag,file.path(DIRT,if(grepl("TRAIN",split_regex))
      "17_MEAN_NESTED_CV_METRICS.csv" else "18_TEST_METRICS.csv"))
    png(file.path(DIRF,if(grepl("TRAIN",split_regex))
      "COMPARACION_MODELOS_TRAIN_NESTED_CV.png" else "COMPARACION_MODELOS_TEST_EXPLORATORIO.png"),
      width=1600,height=1100,res=150)
    ylim <- c(0,1)
    vals <- sapply(FAMILIES,function(fam) {
      z <- ag[ag$family==fam,,drop=FALSE]
      setNames(z$balanced_accuracy,z$panel)[RUN_PANELS]
    })
    palette <- c(CORE="#4E79A7",SHARED="#59A14F",SPECIFIC="#F28E2B",
                 ALL="#B0B0B0",CORE18_47="#9171BD")
    bp <- barplot(vals,beside=TRUE,ylim=ylim,
      col=unname(palette[RUN_PANELS]),
      names.arg=FAMILIES,ylab="Balanced accuracy",xlab="Familia",
      main=if(grepl("TRAIN",split_regex)) "CV anidada condicional en TRAIN" else
        "TEST exploratorio reutilizado — no validacion independiente")
    legend("topright",legend=RUN_PANELS,fill=unname(palette[RUN_PANELS]),bty="n")
    dev.off()
  }
}
# Panel 47 vs CORE del 28: recuento por familia. Si solo se ejecutaron los tres
# nuevos ajustes, las metricas de los 12 anteriores se leen de su corrida previa,
# no se recalculan ni se inventan en esta carpeta.
write_csv(do.call(rbind,lapply(FAMILIES,function(fam) data.frame(
  family=fam, panel28_core=sum(final$family==fam & final$n_questions==3L),
  panel18_core47=nrow(core18_by_family[[fam]]),
  n_overlap_with_panel28=sum(core18_by_family[[fam]]$in_final28),
  n_not_final28=sum(!core18_by_family[[fam]]$in_final28)
))),file.path(DIRT,"21_CORE18_47_VS_CORE28_COUNTS.csv"))
writeLines(c(capture.output(sessionInfo()),"",paste("Fuente 28:",SRC),
             paste("Fuente 18:",consensus18_source)),
  file.path(OUT,"SESSION_INFO.txt"))
cat("\nSCRIPT30_OUTPUT=",OUT,"\n",sep="")
cat("SCRIPT30_STATUS=",paste(bind_or_empty(status)$status,collapse=","),"\n",sep="")

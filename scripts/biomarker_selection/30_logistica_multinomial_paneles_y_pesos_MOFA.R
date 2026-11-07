#!/usr/bin/env Rscript
# METATEST2 / script 30 v1.4.1 — 3 paneles por familia, sin subpaneles.
# MOFA_OTSU (pesos nativos), CORE18_47 y FINAL28 (union final del 28).
# Reutiliza los 6 ajustes ALL/CORE18_47 previos, y ajusta solo 3 MOFA_OTSU.
# No reentrena MOFA, no recalcula Q5/PLSDA/ElasticNet/SGL ni reejecuta el 28.
# Datos: exportacion COMPLETA del script 28 y ready_for_modeling.rds existentes.
# Todos los calculos predictivos son nuevos y condicionados a la seleccion previa.
# Validacion: CV repetida no anidada en TRAIN; TEST descriptivo/reutilizado.
# OVR ROC-AUC, PR average precision por clase; no se usan para seleccionar lambda.

options(stringsAsFactors = FALSE)
CLASS_LEVELS <- c("NP", "FD", "SO")
FAMILIES <- c("P1", "P2", "P3")
VIEWS <- c("tx", "pr", "me", "cl")
PANEL_NAMES <- c("MOFA_OTSU", "CORE18_47", "FINAL28")
RUN_PANELS <- PANEL_NAMES
# Por defecto REUTILIZA los 6 ajustes del ZIP anterior, si validan IDs y RDS.
# MT2_30_REFIT_ALL=1 es la unica forma de recalcular los 9 ajustes.
REFIT_ALL <- identical(Sys.getenv("MT2_30_REFIT_ALL", "0"), "1")
N_CV_REPEATS <- as.integer(Sys.getenv("MT2_30_CV_REPEATS", "5"))
N_OUTER_REPEATS <- as.integer(Sys.getenv("MT2_30_OUTER_REPEATS", "3"))
SEED <- 20261101L
# Smoke: un unico P1/CORE con 1 repeticion de CV estratificada,
# sin validacion anidada imposible por n_FD=3 en TRAIN.
SMOKE <- identical(Sys.getenv("MT2_30_SMOKE", "0"), "1")
# Mantener FAMILIES como el universo de referencia P1/P2/P3 para validar
# las tablas compartidas; limitar solamente el bucle predictivo del smoke.
EXEC_FAMILIES <- if (SMOKE) "P1" else FAMILIES
if (SMOKE) {
  RUN_PANELS <- "MOFA_OTSU"
  N_CV_REPEATS <- 2L
  N_OUTER_REPEATS <- 1L
}
stopifnot(N_CV_REPEATS >= 2, N_OUTER_REPEATS >= 1)

# El script vive en <repo>/scripts/biomarker_selection y detecta la raiz de Git
# a partir de su propia ubicacion (independiente del directorio de trabajo).
SCRIPT30_VERSION <- "1.4.1-3panels-reuse6-new3-feature-importance-comparative-ovr-roc"
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
  "MOFA_OTSU: pesos nativos filtrados por Otsu sobre norma L2 al cuadrado, por familia/vista",
  "CORE18_47: 47 del script 18; FINAL28: todos los biomarcadores finales de cada familia del script 28",
  "El CORE de 33 del script 28 solo es anotacion, NO se ajusta como modelo",
  paste("Paneles ejecutados:",paste(RUN_PANELS,collapse=", ")),
  paste("SMOKE:",SMOKE),
  paste("REFIT_ALL:",REFIT_ALL),
  "Ridge multinomial (alpha=0), regularizacion CV TRAIN; test reutilizado exploratorio",
  "TRAIN: CV estratificada repetida 4-fold; no anidada, porque FD=3 impide CV anidada multinomial glmnet",
  "En tabla 12, n_outer_repeats es el numero de repeticiones OOF (nombre legado, sin folds externos)",
  "Lambda elegida con TRAIN completo (CV 4-fold repetida); OOF en TRAIN reutiliza esa lambda: estimacion optimista, no validacion anidada",
  "TEST reutilizado: exploratorio; ninguno de los metodos constituye validacion externa nueva",
  "ROC y PR: one-vs-rest por clase; AUC ROC y average precision (AP), NO AUC multiclase nativa",
  "La seleccion biologica previa al CV condiciona la CV; no validacion anidada de todo el pipeline",
  "Los betas estandarizados del ranking se calculan como beta_original * SD(TRAIN); importancia L2 entre 3 clases"
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
  ff$feature_label <- if("feature_label" %in% names(e))
     as.character(e$feature_label[ix]) else ff$feature_original
  ff$feature_label[is.na(ff$feature_label) | !nzchar(ff$feature_label)] <-
     ff$feature_original[is.na(ff$feature_label) | !nzchar(ff$feature_label)]
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

# Panel nativo MOFA por familia: todos los predictores por encima de Otsu,
# no solo los que sobrevivieron a Evidence-A/28. No se impone top-N.
MOFA_PANEL <- list(); MOFA_PANEL_AUDIT <- list()
for (fam in FAMILIES) {
  z <- mo[mo$family == fam & !is.na(mo$above_mofa_threshold) &
            mo$above_mofa_threshold, , drop=FALSE]
  assert(nrow(z)>0L, "Sin candidatos MOFA/Otsu para ",fam)
  z$predictor_key <- paste(toupper(z$view_code),z$feature_original,sep="::")
  assert(!anyDuplicated(z$predictor_key),"Claves MOFA repetidas: ",fam)
  # Etiquetas biológicas recuperadas cuando existen en el panel 28 o Evidence-A.
  e <- evidence[evidence$family==fam,,drop=FALSE]
  ekey <- paste(e$view_code,e$feature_original,sep="\r")
  ix <- match(paste(z$view_code,z$feature_original,sep="\r"),ekey)
  z$feature_label <- if("feature_label" %in% names(e)) e$feature_label[ix] else NA_character_
  ix_final <- match(z$key,final$key)
  if("feature_label" %in% names(final)) {
    have_final <- !is.na(ix_final)
    z$feature_label[have_final] <- final$feature_label[ix_final[have_final]]
  }
  z$feature_label[is.na(z$feature_label) | !nzchar(z$feature_label)] <-
    z$feature_original[is.na(z$feature_label) | !nzchar(z$feature_label)]
  MOFA_PANEL[[fam]] <- z
  MOFA_PANEL_AUDIT[[fam]] <- z[,c("family","view_code","feature_original",
                   "predictor_key","MOFA_loading_norm","above_mofa_threshold",
                   "in_final","panel_class","feature_label")]
}
write_csv(do.call(rbind,MOFA_PANEL_AUDIT),file.path(DIRT,"08_MOFA_OTSU_PANEL_MEMBERS.csv"))

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

# Seleccion lambda.1se mediante CV 4-fold repetida sobre TRAIN COMPLETO.
# TRAIN tiene 3 FD: en cada ajuste excluyendo 1 fold quedan al menos 2 FD.
# Una CV anidada (outer+inner) NO es viable porque outer deja solo 2 FD,
# y la inner produciria una unica FD, prohibido por glmnet::lognet.
# El OOF de TRAIN reutiliza una lambda estimada usando TRAIN completo:
# es una evaluacion descriptiva con sesgo optimista (no validacion anidada).
fit_ridge <- function(x,y,seed,nrep=N_CV_REPEATS) {
  y <- factor(y, levels=CLASS_LEVELS)
  assert(all(table(y)>=2L), "El ajuste requiere al menos 2 observaciones por clase: ",
         paste(names(table(y)),as.integer(table(y)),collapse=", "))
  k <- 4L
  assert(nrow(x)>=k && min(table(y))>=2L,
         "TRAIN insuficiente para CV interna estratificada de 4 folds")
  cv <- lapply(seq_len(nrep),function(i) {
    foldid <- stratified_folds(y,k,seed+i)
    assert(length(unique(foldid))==k && all(tabulate(foldid,nbins=k)>0L),
           "Folds CV internos incompletos")
    for(j in seq_len(k)) {
      train_fold <- y[foldid!=j]
      assert(all(table(train_fold)>0),
             "Una clase ausente en entrenamiento de CV interna")
    }
    glmnet::cv.glmnet(x=x,y=y,family="multinomial",
      type.multinomial="grouped",alpha=0,standardize=TRUE,
      weights=class_weights(y), foldid=foldid, nfolds=k,
      type.measure="deviance")
  })
  lambda_1se <- vapply(cv,function(z)as.numeric(z$lambda.1se),numeric(1))
  assert(all(is.finite(lambda_1se) & lambda_1se>0),
         "cv.glmnet no devolvio lambda.1se valida")
  lambda <- stats::median(lambda_1se)
  fit <- glmnet::glmnet(x=x,y=y,family="multinomial",
      type.multinomial="grouped",alpha=0,standardize=TRUE,
      weights=class_weights(y))
  list(fit=fit,lambda=lambda,lambda_1se=lambda_1se,
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
# ROC y PR one-vs-rest calculadas a partir de probabilidades de prediccion.
# Puntos agrupados por score para tratar empates sin introducir orden artificial.
# PR: average precision (AP) escalonada; no se denomina ROC-PR AUC indistintamente.
# Si falta una clase positiva o la clase negativa, se informa NA, no 0.
binary_ovr <- function(truth_binary,score,fam,panel,split,cls) {
  keep <- is.finite(score) & !is.na(truth_binary)
  y <- as.logical(truth_binary[keep]); s <- as.numeric(score[keep])
  pos <- sum(y); neg <- sum(!y); n <- length(y)
  metric <- data.frame(family=fam,panel=panel,split=split,class=cls,
    n=n,n_pos=pos,n_neg=neg,prevalence=if(n>0L)pos/n else NA_real_,
    roc_auc=NA_real_,average_precision=NA_real_,stringsAsFactors=FALSE)
  if(pos==0L || neg==0L) return(list(metric=metric,roc=NULL,pr=NULL))
  ord <- order(-s,seq_along(s))
  y <- y[ord]; s <- s[ord]
  ends <- c(which(diff(s)!=0L),length(s))
  tp <- c(0L,cumsum(y)[ends]); fp <- c(0L,ends-cumsum(y)[ends])
  tpr <- tp/pos; fpr <- fp/neg
  roc_auc <- sum(diff(fpr)*(head(tpr,-1L)+tail(tpr,-1L))/2)
  precision <- c(1,tp[-1L]/(tp[-1L]+fp[-1L]))
  recall <- tpr
  ap <- sum(diff(recall)*precision[-1L])
  metric$roc_auc <- roc_auc; metric$average_precision <- ap
  roc <- data.frame(family=fam,panel=panel,split=split,class=cls,
    fpr=fpr,tpr=tpr,threshold=c(Inf,s[ends]),stringsAsFactors=FALSE)
  pr <- data.frame(family=fam,panel=panel,split=split,class=cls,
    recall=recall,precision=precision,threshold=c(Inf,s[ends]),
    prevalence=pos/n,stringsAsFactors=FALSE)
  list(metric=metric,roc=roc,pr=pr)
}
discrimination_ovr <- function(probs,y,fam,panel,split) {
  assert(all(is.finite(probs)),"Probabilidades no finitas: ",fam,"/",panel,"/",split)
  assert(all(CLASS_LEVELS %in% colnames(probs)),"Faltan columnas probabilidad")
  assert(all(abs(rowSums(probs)-1)<1e-6),"Probabilidades sin suma 1")
  z <- lapply(CLASS_LEVELS,function(cls)
    binary_ovr(as.character(y)==cls,probs[,cls],fam,panel,split,cls))
  rows <- do.call(rbind,lapply(z,`[[`,"metric"))
  join <- function(field) {
    a <- lapply(z,`[[`,field)
    a <- Filter(Negate(is.null),a)
    if(length(a))do.call(rbind,a) else NULL
  }
  list(metrics=rows,roc=join("roc"),pr=join("pr"))
}
plot_ovr_curve <- function(curves,summary,kind,fam,panel,split,destination) {
  if(is.null(curves) || !nrow(curves)) return(invisible(FALSE))
  isroc <- identical(kind,"ROC")
  grDevices::png(destination,width=1650,height=1300,res=165,bg="white")
  on.exit(grDevices::dev.off(),add=TRUE)
  graphics::par(mar=c(5.2,5.5,4.5,1.5))
  graphics::plot(NA_real_,xlim=c(0,1),ylim=c(0,1),xlab=if(isroc)"False positive rate" else "Recall",
    ylab=if(isroc)"True positive rate" else "Precision",
    main=paste(fam,panel,split,if(isroc)"ROC one-vs-rest" else "Precision-Recall one-vs-rest",sep=" | "),
    xaxs="i",yaxs="i")
  if(isroc) graphics::abline(a=0,b=1,lty=3,col="grey55")
  cc <- c(NP="#4E79A7",FD="#E15759",SO="#59A14F")
  legends <- character(); colors <- character()
  for(cls in CLASS_LEVELS) {
    z <- curves[curves$class==cls,,drop=FALSE]
    if(!nrow(z))next
    o <- summary[summary$class==cls,,drop=FALSE]
    if(isroc) graphics::lines(z$fpr,z$tpr,col=cc[cls],lwd=3,type="s")
    else {
      graphics::lines(z$recall,z$precision,col=cc[cls],lwd=3,type="s")
      graphics::abline(h=o$prevalence[1],col=cc[cls],lty=3)
    }
    score <- if(isroc)o$roc_auc[1] else o$average_precision[1]
    legends <- c(legends,sprintf("%s (n+=%d; %s=%.3f)",cls,o$n_pos[1],
      if(isroc)"AUC" else "AP",score))
    colors <- c(colors,cc[cls])
  }
  if(length(legends)) graphics::legend("bottomright",legend=legends,col=colors,
    lwd=3,bty="n",cex=0.83)
  graphics::mtext(if(grepl("^TEST",split))
    "TEST reutilizado: exploratorio; FD tiene muy pocas observaciones" else
    "TRAIN CV no anidada: lambda y biomarcadores seleccionados con TRAIN",
    side=1,line=3.7,cex=0.72)
  invisible(TRUE)
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
  discr <- discrimination_ovr(probs,y,fam,panel,split)
  list(metrics=metr,classes=classes,pred=pred,conf=conf,discr=discr)
}

metrics <- list(); class_metrics <- list(); predictions <- list(); confusions <- list()
coefficients <- list(); eligibility <- list(); configs <- list(); status <- list()
discrimination <- list(); roc_points <- list(); pr_points <- list()
register_discrimination <- function(id,res,fam,panel,split,make_plot=FALSE) {
  discr <- res$discr
  discrimination[[id]] <<- discr$metrics
  roc_points[[id]] <<- discr$roc
  pr_points[[id]] <<- discr$pr
  if(make_plot) {
    for(kind in c("ROC","PR")) {
      points <- if(kind=="ROC")discr$roc else discr$pr
      plot_ovr_curve(points,discr$metrics,kind,fam,panel,split,
        file.path(DIRF,paste0(kind,"_OVR_",fam,"_",panel,"_",split,".png")))
    }
  }
}


# Reutilizar los 6 modelos ya ajustados por la v1.3: ALL -> FINAL28,
# CORE18_47 -> CORE18_47. Comprobar artefactos, IDs, rutas y lambda antes.
find_prior_successful <- function() {
  root <- dirname(OUT)
  dirs <- list.dirs(root, recursive=FALSE,full.names=TRUE)
  dirs <- setdiff(dirs, OUT)
  dirs <- dirs[order(file.info(dirs)$mtime,decreasing=TRUE)]
  for (d in dirs) {
    stp <- file.path(d,"tables","10_MODEL_STATUS.csv")
    cfgp <- file.path(d,"tables","12_MODEL_HYPERPARAMETERS.csv")
    if(!file.exists(stp) || !file.exists(cfgp)) next
    st <- tryCatch(read_csv(stp,c("family","panel","status")),error=function(e)NULL)
    if(is.null(st)) next
    need <- expand.grid(family=FAMILIES,panel=c("ALL","CORE18_47"),stringsAsFactors=FALSE)
    ok <- vapply(seq_len(nrow(need)),function(i) {
      q <- st[st$family==need$family[i] & st$panel==need$panel[i],,drop=FALSE]
      nrow(q)==1L && identical(q$status,"PASS") &&
        file.exists(file.path(d,"models",paste0(need$family[i],"_",need$panel[i],"_ridge_multinomial.rds")))
    },logical(1))
    if(all(ok)) return(d)
  }
  NA_character_
}
PRIOR_RUN <- if(REFIT_ALL || SMOKE) NA_character_ else find_prior_successful()
if (!REFIT_ALL && !SMOKE) {
  assert(!is.na(PRIOR_RUN),
    "No se encontro corrida anterior con ALL y CORE18_47 PASS en las 3 familias. ",
    "No se recalcularan sin permiso. Para forzar modelos nuevos: MT2_30_REFIT_ALL=1")
  message("SEIS_MODELOS_REUTILIZADOS_DESDE=",PRIOR_RUN)
  old_readme <- readLines(file.path(PRIOR_RUN,"README.txt"),warn=FALSE)
  old_source <- sub("^Script 28: ","",grep("^Script 28: ",old_readme,value=TRUE))
  assert(length(old_source)==1L &&
         identical(normalizePath(old_source,mustWork=TRUE),normalizePath(SRC,mustWork=TRUE)),
         "La fuente de 28 difiere del ajuste anterior; se prohibe reutilizar")
  old_18 <- readLines(file.path(PRIOR_RUN,"CONSENSO_47_ORIGEN.txt"),warn=FALSE)
  old_consensus <- sub("^Consenso del script 18: ","",
    grep("^Consenso del script 18: ",old_18,value=TRUE))
  assert(length(old_consensus)==1L &&
         identical(normalizePath(old_consensus,mustWork=TRUE),
                   normalizePath(consensus18_source,mustWork=TRUE)),
         "El panel de 47 difiere de la referencia anterior; se prohibe reutilizar")
  writeLines(c(paste("PRIOR_RUN:",normalizePath(PRIOR_RUN)),
    "Se reutilizan exactamente ALL -> FINAL28 y CORE18_47 para P1/P2/P3",
    "Se leen RDS, coeficientes, matrices de confusion, predicciones y curvas preexistentes."),
    file.path(OUT,"ORIGEN_MODELOS_REUTILIZADOS.txt"))
}
read_prior_panel <- function(fname, fam, old_panel, new_panel, required=character()) {
  p <- file.path(PRIOR_RUN,"tables",fname)
  x <- read_csv(p,c("family","panel",required))
  z <- x[x$family==fam & x$panel==old_panel,,drop=FALSE]
  if(nrow(z)) z$panel <- new_panel
  z
}
# Copiar tablas de evaluacion/prediccion de los 6 modelos ya ajustados.
reuse_model <- function(fam,panel,dd,ff,key) {
  old_panel <- if(panel=="FINAL28") "ALL" else panel
  prior_rds <- file.path(PRIOR_RUN,"models",
                          paste0(fam,"_",old_panel,"_ridge_multinomial.rds"))
  mod <- readRDS(prior_rds)
  assert(identical(mod$predictor_keys,colnames(dd$xtr)),
         "Predictores/reordenacion diferentes al RDS anterior: ",fam,"/",panel)
  assert(identical(as.integer(mod$n_train),as.integer(nrow(dd$xtr))) &&
         identical(as.character(mod$family),fam) &&
         identical(as.character(mod$panel),old_panel),
         "Datos/familia diferentes al modelo previo: ",key)
  prior_cfg <- read_prior_panel("12_MODEL_HYPERPARAMETERS.csv",fam,old_panel,panel,
                               c("n_selected","n_eligible","lambda","ready_for_modeling"))
  assert(nrow(prior_cfg)==1L && prior_cfg$n_selected==nrow(ff) &&
         prior_cfg$n_eligible==ncol(dd$xtr) &&
         identical(normalizePath(prior_cfg$ready_for_modeling,mustWork=TRUE),
                   normalizePath(dat$rds_path,mustWork=TRUE)) &&
         isTRUE(all.equal(as.numeric(prior_cfg$lambda),as.numeric(mod$lambda))),
         "La referencia previa no coincide con el panel actual: ",key)
  out_model <- file.path(OUT,"models",paste0(key,"_ridge_multinomial.rds"))
  assert(file.copy(prior_rds,out_model,overwrite=FALSE),
         "No se pudo copiar RDS reutilizado: ",key)
  configs[[key]] <<- prior_cfg
  # Recuperar sin reentrenar TODOS los puntos ROC/PR y las predicciones.
  files <- list(
    metrics=list("13_MODEL_METRICS.csv",c("split","accuracy")),
    class_metrics=list("14_PER_CLASS_METRICS.csv",c("split","class")),
    predictions=list("15_PREDICTIONS_PER_SAMPLE.csv",c("split","sample_id")),
    coefficients=list("16_LOGISTIC_COEFFICIENTS.csv",c("class","predictor_key")),
    confusions=list("19_CONFUSION_LONG.csv",c("split")),
    discrimination=list("25_ROC_PR_OVR_BY_CLASS.csv",c("split","class")),
    roc_points=list("23_ROC_CURVE_POINTS_OVR.csv",c("split","class")),
    pr_points=list("24_PR_CURVE_POINTS_OVR.csv",c("split","class")))
  for(nm in names(files)) {
    tbl <- read_prior_panel(files[[nm]][[1]],fam,old_panel,panel,files[[nm]][[2]])
    assert(nrow(tbl)>0L,"Tabla reutilizada vacia: ",nm," ",key)
    # asignar al contenedor global adecuado, sin confundir con variables locales.
    assign(nm, c(get(nm,envir=.GlobalEnv),setNames(list(tbl),key)),envir=.GlobalEnv)
  }
  # Regenerar figuras a partir de PUNTOS GUARDADOS: no recalcular probabilidades.
  rr <- roc_points[[key]]; pp <- pr_points[[key]]; summary <- discrimination[[key]]
  splits <- unique(as.character(summary$split))
  splits <- splits[grepl("^TEST",splits) | grepl("REPEAT_1$",splits)]
  for(sp in splits) {
    z <- summary[summary$split==sp,,drop=FALSE]
    for(kind in c("ROC","PR")) {
      points <- if(kind=="ROC")rr else pp
      points <- points[points$split==sp,,drop=FALSE]
      plot_ovr_curve(points,z,kind,fam,panel,sp,
        file.path(DIRF,paste0(kind,"_OVR_",fam,"_",panel,"_",sp,".png")))
    }
  }
  status[[key]] <<- data.frame(family=fam,panel=panel,status="PASS_REUSED",
    detail=paste0("RDS+predicciones+ROC/PR reutilizados desde ",basename(PRIOR_RUN)))
  message("REUSADO: ",key," (lambda=",signif(mod$lambda,4),")")
  invisible(mod)
}

importance_rows <- list(); importance_class_rows <- list()
plot_importance <- function(z, family, panel) {
  if(!nrow(z)) return(invisible(NULL))
  # Solo el grafico limita el numero de etiquetas; CSV incluye TODAS las features.
  zp <- head(z[order(-z$importance_l2_standardized,z$predictor_key),,drop=FALSE],30L)
  zp <- zp[order(zp$importance_l2_standardized),,drop=FALSE]
  png(file.path(DIRF,paste0("IMPORTANCE_",family,"_",panel,"_TOP30.png")),
      width=2200,height=max(1000,85*nrow(zp)),res=165,bg="white")
  par(mar=c(6,11,5,3))
  barplot(zp$importance_l2_standardized,horiz=TRUE,las=1,
          names.arg=zp$feature_label,cex.names=0.72,border=NA,
          col=c(tx="#B07AA1",pr="#4E79A7",me="#59A14F",cl="#E15759")[zp$view_code],
          xlab="Norma L2 de coeficientes ridge estandarizados (NP/FD/SO)",
          main=paste(family,panel,"- importancia regresion"))
  legend("bottomright",legend=c("TX","PR","ME","CL"),
         fill=c("#B07AA1","#4E79A7","#59A14F","#E15759"),bty="n")
  dev.off()
  for(cls in CLASS_LEVELS) {
    cc <- paste0("beta_std_",cls)
    zz <- head(z[order(-abs(z[[cc]]),z$feature_label),,drop=FALSE],20L)
    zz <- zz[order(zz[[cc]]),,drop=FALSE]
    png(file.path(DIRF,paste0("BETA_",family,"_",panel,"_",cls,"_TOP20.png")),
        width=2000,height=1500,res=160,bg="white")
    par(mar=c(5,11,5,2))
    barplot(zz[[cc]],horiz=TRUE,las=1,names.arg=zz$feature_label,
      col=ifelse(zz[[cc]]>=0,"#2370B7","#D55E00"),border=NA,cex.names=0.7,
      xlab="Coeficiente estandarizado de score multinomial (signo orientativo)",
      main=paste(family,panel,cls,"- top 20 por |beta|"))
    abline(v=0,col="black")
    dev.off()
  }
}
for(fam in EXEC_FAMILIES) {
  message("Leyendo matrices ya preprocesadas: ",fam)
  dat <- load_family(fam,bio_run)
  for(panel in RUN_PANELS) {
    ff <- switch(panel,
      MOFA_OTSU=MOFA_PANEL[[fam]],
      CORE18_47=core18_by_family[[fam]],
      FINAL28=final[final$family==fam,,drop=FALSE])
    key <- paste(fam,panel,sep="_")
    message("Modelo ",key," (",nrow(ff)," features de panel)")
    if(!nrow(ff)) {
      status[[key]] <- data.frame(family=fam,panel=panel,
                                 status="SKIP_EMPTY_PANEL",detail="Sin features")
      next
    }
    dd <- matrices_for_panel(dat,ff,fam,panel)
    eligibility[[key]] <- dd$eligibility
    if(!is.null(dd$skip)) {
      status[[key]] <- data.frame(family=fam,panel=panel,
                                status="SKIP_NO_FEATURES",detail=dd$skip)
      next
    }
    reused <- !REFIT_ALL && !SMOKE && panel %in% c("CORE18_47","FINAL28")
    if(reused) {
      reuse_model(fam,panel,dd,ff,key)
    } else {
      fitted <- tryCatch(fit_ridge(dd$xtr,dd$ytr,
             SEED+match(fam,FAMILIES)*100+match(panel,PANEL_NAMES)),
             error=function(e)e)
      if(inherits(fitted,"error")) {
        status[[key]] <- data.frame(family=fam,panel=panel,
          status="FAIL_FINAL_FIT",detail=conditionMessage(fitted)); next
      }
      err <- NULL
      for(rep in seq_len(N_OUTER_REPEATS)) {
        foldid <- stratified_folds(dd$ytr,4L,SEED+10000L*match(fam,FAMILIES)+
             1000L*match(panel,PANEL_NAMES)+rep)
        assert(all(tabulate(foldid,nbins=4L)>0L),"Fold OOF vacio")
        pred_rep <- matrix(NA_real_,nrow(dd$xtr),3L,
          dimnames=list(rownames(dd$xtr),CLASS_LEVELS))
        for(k in 1:4) {
          itr <- which(foldid!=k); ite <- which(foldid==k)
          tryCatch({
            yy <- droplevels(dd$ytr[itr])
            assert(all(table(yy)>=2L),"Insuficientes observaciones por clase en OOF")
            f_oo <- glmnet::glmnet(x=dd$xtr[itr,,drop=FALSE],y=yy,
              family="multinomial",type.multinomial="grouped",alpha=0,
              standardize=TRUE,weights=class_weights(yy),lambda=fitted$lambda)
            pred_rep[ite,] <- prob_matrix(f_oo,dd$xtr[ite,,drop=FALSE],fitted$lambda)
          },error=function(e){err <<- paste0("rep=",rep," fold=",k,": ",conditionMessage(e))})
          if(!is.null(err))break
        }
        if(!is.null(err))break
        assert(all(is.finite(pred_rep)),"Predicciones no finitas CV TRAIN")
        sn <- paste0("TRAIN_REPEATED_CV_TUNED_FULL_TRAIN_REPEAT_",rep)
        res <- assess(pred_rep,dd$ytr,fam,panel,sn)
        metrics[[paste0(key,"_CV_",rep)]] <- res$metrics
        class_metrics[[paste0(key,"_CV_",rep)]] <- res$classes
        predictions[[paste0(key,"_CV_",rep)]] <- res$pred
        confusions[[paste0(key,"_CV_",rep)]] <- data.frame(family=fam,panel=panel,
          split=sn,as.data.frame(res$conf))
        register_discrimination(paste0(key,"_CV_",rep),res,fam,panel,sn,
          make_plot=(rep==1L))
      }
      if(!is.null(err)) {
        status[[key]] <- data.frame(family=fam,panel=panel,
          status="FAIL_TRAIN_REPEATED_CV",detail=err); next
      }
      test_p <- prob_matrix(fitted$fit,dd$xte,fitted$lambda)
      res <- assess(test_p,dd$yte,fam,panel,"TEST_REUSED_EXPLORATORY")
      metrics[[paste0(key,"_TEST")]] <- res$metrics
      class_metrics[[paste0(key,"_TEST")]] <- res$classes
      predictions[[paste0(key,"_TEST")]] <- res$pred
      confusions[[paste0(key,"_TEST")]] <- data.frame(family=fam,panel=panel,
        split="TEST_REUSED_EXPLORATORY",as.data.frame(res$conf))
      register_discrimination(paste0(key,"_TEST"),res,fam,panel,
                              "TEST_REUSED_EXPLORATORY",make_plot=TRUE)
      saveRDS(list(model=fitted$fit,lambda=fitted$lambda,
        predictor_keys=colnames(dd$xtr),family=fam,panel=panel,
        n_train=nrow(dd$xtr),classes=CLASS_LEVELS),
        file.path(OUT,"models",paste0(key,"_ridge_multinomial.rds")))
      cf <- stats::coef(fitted$fit,s=fitted$lambda)
      coefficients[[key]] <- do.call(rbind,lapply(names(cf),function(cls) {
        cc <- as.matrix(cf[[cls]])
        data.frame(family=fam,panel=panel,class=cls,
          predictor_key=rownames(cc),coefficient=as.numeric(cc[,1L]))
      }))
      configs[[key]] <- data.frame(family=fam,panel=panel,n_train=nrow(dd$xtr),
        n_test=nrow(dd$xte),n_selected=nrow(ff),n_eligible=ncol(dd$xtr),
        lambda=fitted$lambda,lambda_1se_median=median(fitted$lambda_1se),
        cv_deviance_median=median(fitted$cv_deviance),alpha=0,
        n_cv_repeats=N_CV_REPEATS,n_outer_repeats=N_OUTER_REPEATS,
        ready_for_modeling=dat$rds_path)
      status[[key]] <- data.frame(family=fam,panel=panel,
          status="PASS",detail="NEW_FIT: FINAL TRAIN + TEST reutilizado")
    }
    # Importancia de cada predictor: beta(glmnet,escala original) * SD(TRAIN).
    # Norma de beta estandarizado para comparar dentro de UN modelo/familia.
    cf <- coefficients[[key]]
    assert(nrow(cf)>0L,"Sin coeficientes en ",key)
    ii <- seq_along(colnames(dd$xtr)); keys <- colnames(dd$xtr)
    sdtr <- apply(dd$xtr,2L,stats::sd)
    bet <- sapply(CLASS_LEVELS,function(cls) {
      z <- cf[cf$class==cls,,drop=FALSE]
      assert(!anyDuplicated(z$predictor_key),"Coeficientes duplicados: ",key,cls)
      ix <- match(keys,z$predictor_key)
      assert(!anyNA(ix),"Falta beta para ",cls," en ",key)
      as.numeric(z$coefficient[ix])
    })
    if(is.null(dim(bet))) bet <- matrix(bet,ncol=length(CLASS_LEVELS))
    colnames(bet) <- CLASS_LEVELS
    bstd <- sweep(bet,1L,sdtr,"*")
    ix <- match(keys,ff$predictor_key)
    assert(!anyNA(ix),"Feature sin metadatos de panel: ",key)
    z <- data.frame(family=fam,panel=panel,predictor_key=keys,
      view_code=ff$view_code[ix],feature_original=ff$feature_original[ix],
      feature_label=if("feature_label" %in% names(ff))
                    ff$feature_label[ix] else ff$feature_original[ix],
      sd_train=as.numeric(sdtr),beta_std_NP=bstd[,"NP"],
      beta_std_FD=bstd[,"FD"],beta_std_SO=bstd[,"SO"],
      importance_l2_standardized=sqrt(rowSums(bstd^2)),
      stringsAsFactors=FALSE)
    lookup <- match(paste(fam,z$view_code,z$feature_original,sep="\r"),mo$key)
    z$MOFA_loading_norm <- mo$MOFA_loading_norm[lookup]
    z$MOFA_above_otsu <- mo$above_mofa_threshold[lookup]
    z$in_final28 <- mo$in_final[lookup]
    z <- z[order(-z$importance_l2_standardized,z$predictor_key),,drop=FALSE]
    z$rank_overall <- seq_len(nrow(z))
    z$rank_within_view <- ave(z$importance_l2_standardized,z$view_code,
                              FUN=function(x)rank(-x,ties.method="min"))
    importance_rows[[key]] <- z
    importance_class_rows[[key]] <- do.call(rbind,lapply(CLASS_LEVELS,function(cls){
      b <- z[[paste0("beta_std_",cls)]]
      data.frame(family=fam,panel=panel,class=cls,predictor_key=z$predictor_key,
        view_code=z$view_code,feature_label=z$feature_label,
        beta_standardized=b,abs_beta_standardized=abs(b),
        rank=rank(-abs(b),ties.method="min"))
    }))
    plot_importance(z,fam,panel)
  }
}
write_csv(if(length(importance_rows))do.call(rbind,importance_rows) else
  data.frame(),file.path(DIRT,"30_LOGISTIC_FEATURE_IMPORTANCE_GLOBAL.csv"))
write_csv(if(length(importance_class_rows))do.call(rbind,importance_class_rows) else
  data.frame(),file.path(DIRT,"31_LOGISTIC_FEATURE_IMPORTANCE_BY_CLASS.csv"))
# Comparar los rankings nativos MOFA y logísticos; Spearman descriptivo por vista.
if(length(importance_rows)) {
  rank_comp <- do.call(rbind,lapply(importance_rows,function(z) {
    do.call(rbind,lapply(unique(z$view_code),function(v) {
      zz <- z[z$view_code==v & is.finite(z$MOFA_loading_norm),,drop=FALSE]
      data.frame(family=z$family[1],panel=z$panel[1],view_code=v,
        n_features=nrow(zz),
        spearman_mofa_vs_logistic=if(nrow(zz)>=3L &&
          length(unique(zz$MOFA_loading_norm))>1L &&
          length(unique(zz$importance_l2_standardized))>1L)
          cor(zz$MOFA_loading_norm,zz$importance_l2_standardized,
              method="spearman") else NA_real_)
    }))
  }))
  write_csv(rank_comp,file.path(DIRT,"32_MOFA_VS_LOGISTIC_RANK_CORRELATION.csv"))
  panel_audit <- do.call(rbind,lapply(importance_rows,function(z)
    data.frame(family=z$family[1],panel=z$panel[1],
      n_eligible=nrow(z),
      n_above_mofa_threshold=sum(z$MOFA_above_otsu,na.rm=TRUE),
      n_in_final28=sum(z$in_final28,na.rm=TRUE))))
  write_csv(panel_audit,file.path(DIRT,"33_PANEL_SIZES_AND_OVERLAPS.csv"))
}

bind_or_empty <- function(x,cols=character()) {
  if(length(x)) return(do.call(rbind,x))
  as.data.frame(setNames(rep(list(character(0)),length(cols)),cols),
                check.names=FALSE,stringsAsFactors=FALSE)
}
write_csv(bind_or_empty(status,c("family","panel","status","detail")),file.path(DIRT,"10_MODEL_STATUS.csv"))
write_csv(bind_or_empty(eligibility,c("family","panel","predictor_key","view_code","feature_original","reason")),file.path(DIRT,"11_MODEL_FEATURE_ELIGIBILITY.csv"))
write_csv(bind_or_empty(configs,c("family","panel","n_train","n_test","n_selected","n_eligible","lambda","lambda_1se_median","cv_deviance_median","alpha","n_cv_repeats","n_outer_repeats","ready_for_modeling")),file.path(DIRT,"12_MODEL_HYPERPARAMETERS.csv"))
write_csv(bind_or_empty(metrics,c("family","panel","split","n","accuracy","balanced_accuracy","macro_f1","log_loss","brier")),file.path(DIRT,"13_MODEL_METRICS.csv"))
write_csv(bind_or_empty(class_metrics,c("family","panel","split","class","n","recall","precision","f1")),file.path(DIRT,"14_PER_CLASS_METRICS.csv"))
write_csv(bind_or_empty(predictions,c("family","panel","split","sample_id","true","predicted","prob_NP","prob_FD","prob_SO")),file.path(DIRT,"15_PREDICTIONS_PER_SAMPLE.csv"))
write_csv(bind_or_empty(coefficients,c("family","panel","class","predictor_key","coefficient")),file.path(DIRT,"16_LOGISTIC_COEFFICIENTS.csv"))
write_csv(bind_or_empty(confusions,c("family","panel","split","truth","predicted","Freq")),file.path(DIRT,"19_CONFUSION_LONG.csv"))
write_csv(bind_or_empty(discrimination,c("family","panel","split","class","n","n_pos","n_neg","prevalence","roc_auc","average_precision")),file.path(DIRT,"25_ROC_PR_OVR_BY_CLASS.csv"))
write_csv(bind_or_empty(roc_points,c("family","panel","split","class","fpr","tpr","threshold")),file.path(DIRT,"23_ROC_CURVE_POINTS_OVR.csv"))
write_csv(bind_or_empty(pr_points,c("family","panel","split","class","recall","precision","threshold","prevalence")),file.path(DIRT,"24_PR_CURVE_POINTS_OVR.csv"))
if(length(discrimination)) {
  dd <- bind_or_empty(discrimination)
  by_split <- split(dd,interaction(dd$family,dd$panel,dd$split,drop=TRUE))
  discr_summary <- do.call(rbind,lapply(by_split,function(z) data.frame(
    family=z$family[1],panel=z$panel[1],split=z$split[1],
    n_classes_evaluable=sum(is.finite(z$roc_auc)),
    macro_roc_auc=if(any(is.finite(z$roc_auc))) mean(z$roc_auc,na.rm=TRUE) else NA_real_,
    macro_average_precision=if(any(is.finite(z$average_precision)))
      mean(z$average_precision,na.rm=TRUE) else NA_real_)))
  write_csv(discr_summary,file.path(DIRT,"26_ROC_PR_MACRO_BY_MODEL.csv"))
}

if(length(metrics)) {
  mm <- bind_or_empty(metrics)
  for(split_regex in c("^TRAIN_REPEATED_CV","^TEST_REUSED")) {
    ss <- mm[grepl(split_regex,mm$split),,drop=FALSE]
    if(!nrow(ss)) next
    ag <- stats::aggregate(cbind(accuracy,balanced_accuracy,macro_f1) ~ family + panel,
      ss,FUN=mean)
    write_csv(ag,file.path(DIRT,if(grepl("TRAIN",split_regex))
      "17_MEAN_REPEATED_CV_METRICS.csv" else "18_TEST_METRICS.csv"))
    png(file.path(DIRF,if(grepl("TRAIN",split_regex))
      "COMPARACION_MODELOS_TRAIN_REPEATED_CV.png" else "COMPARACION_MODELOS_TEST_EXPLORATORIO.png"),
      width=1600,height=1100,res=150)
    ylim <- c(0,1)
    vals <- vapply(EXEC_FAMILIES,function(fam) {
      z <- ag[ag$family==fam,,drop=FALSE]
      unname(setNames(z$balanced_accuracy,z$panel)[RUN_PANELS])
    },numeric(length(RUN_PANELS)))
    vals <- matrix(vals, nrow=length(RUN_PANELS), ncol=length(EXEC_FAMILIES),
                   dimnames=list(RUN_PANELS, EXEC_FAMILIES))
    palette <- c(MOFA_OTSU="#E15759",CORE18_47="#9171BD",FINAL28="#4E79A7")
    bp <- barplot(vals,beside=TRUE,ylim=ylim,
      col=unname(palette[RUN_PANELS]),
      names.arg=EXEC_FAMILIES,ylab="Balanced accuracy",xlab="Familia",
      main=if(grepl("TRAIN",split_regex)) "CV repetida TRAIN (lambda calibrada en mismos datos)" else
        "TEST exploratorio reutilizado — no validacion independiente")
    legend("topright",legend=RUN_PANELS,fill=unname(palette[RUN_PANELS]),bty="n")
    dev.off()
  }
}
# Panel 47 vs CORE (33) del 28: tabla auxiliar de solapamiento, sin ajuste CORE.
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
status_table <- bind_or_empty(status)
if(!any(grepl("^PASS",status_table$status)) || any(!grepl("^PASS",status_table$status))) {
  stop("SCRIPT30_INCOMPLETE: modelos sin ajuste correcto. Revisar ",
       file.path(DIRT,"10_MODEL_STATUS.csv"),call.=FALSE)
}

# Comparacion ROC TEST por clase: tres curvas de los modelos para cada familia.
# Reutiliza exclusivamente los puntos ROC ya exportados; no ajusta modelos.
if (!SMOKE) {
  roc_comparison_script <- file.path(dirname(script_file),
                                     "30b_ROC_comparativa_modelos_OVR.R")
  assert(file.exists(roc_comparison_script),
         "Falta script de figuras comparativas: ", roc_comparison_script)
  source(roc_comparison_script, local = new.env(parent = globalenv()), echo = FALSE)
}

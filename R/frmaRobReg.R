frmaRobReg <- function(object, background, normalize, summarize, target, input.vecs, output.param, verbose){

  if(class(object)=="AffyBatch"){
    cdfname <- cleancdfname(cdfName(object))
    tmp <- gsub("cdf","",cdfname)
    platform <- gsub("..entrezg", "", tmp)
    vecdataname <- paste(tmp, "frmavecs", sep="")
  }
  if(class(object)%in%c("ExonFeatureSet","GeneFeatureSet")){
    cdfname <- annotation(object)
    platform <- gsub("pd.","",cdfname,fixed=TRUE)
    vecdataname <- paste(platform, "frmavecs", sep="")
  }
  
  if(background == "rma"){
    if(verbose) message("Background Correcting ...\n")
    if(class(object)=="AffyBatch") object <- bg.correct.rma(object)
    if(class(object)%in%c("ExonFeatureSet","GeneFeatureSet")) object <- backgroundCorrect(object, verbose=FALSE)
  }

  if(class(object)%in%c("ExonFeatureSet","GeneFeatureSet")){
    featureInfo <- getFidProbeset(object)
    pmi <- featureInfo[["fid"]]
    pns <- as.character(featureInfo[["fsetid"]])
    pms <- exprs(object)[pmi,, drop=FALSE]
  }
  if(class(object)=="AffyBatch"){
    pms <- pm(object)
    pns <- probeNames(object)
  }

  if( is.null(input.vecs$normVec) | is.null(input.vecs$probeVec) |
     is.null(input.vecs$probeVarWithin) | is.null(input.vecs$probeVarBetween) |
     (summarize=="robust_weighted_average" & is.null(input.vecs$probesetSD)) |
     ((is.null(input.vecs$probeVecCore)|is.null(input.vecs$mapCore)) & target=="core") |
     ((is.null(input.vecs$probeVecExt)|is.null(input.vecs$mapExt)) & target=="extended") |
     ((is.null(input.vecs$probeVecFull)|is.null(input.vecs$mapFull)) & target=="full") ){
    pkg <- paste(platform, "frmavecs", sep="")
    require(pkg, character.only=TRUE, quietly=TRUE) || stop(paste(pkg, "package must be installed first"))
    data(list=eval(vecdataname))

    if(is.null(input.vecs$normVec)) input.vecs$normVec <- get(vecdataname)$normVec
    if(is.null(input.vecs$probeVec)) input.vecs$probeVec <- get(vecdataname)$probeVec
    if(is.null(input.vecs$probeVarWithin)) input.vecs$probeVarWithin <- get(vecdataname)$probeVarWithin
    if(is.null(input.vecs$probeVarBetween)) input.vecs$probeVarBetween <- get(vecdataname)$probeVarBetween
    if((is.null(input.vecs$probeVecCore)|is.null(input.vecs$mapCore)) & target=="core"){
      input.vecs$probeVecCore <- get(vecdataname)$probeVecCore
      input.vecs$mapCore <- get(vecdataname)$mapCore
    }
    if((is.null(input.vecs$probeVecExt)|is.null(input.vecs$mapExt)) & target=="extended"){
      input.vecs$probeVecExt <- get(vecdataname)$probeVecExt
      input.vecs$mapExt <- get(vecdataname)$mapExt
    }
    if((is.null(input.vecs$probeVecFull)|is.null(input.vecs$mapFull)) & target=="full"){
      input.vecs$probeVecFull <- get(vecdataname)$probeVecFull
      input.vecs$mapFull <- get(vecdataname)$mapFull
    }
    if(is.null(input.vecs$probesetSD) & summarize=="robust_weighted_average" & (class(object)=="AffyBatch" | target=="probeset")) input.vecs$probesetSD <- get(vecdataname)$probesetSD
  }

  if(normalize == "quantile"){
    if(verbose) message("Normalizing ...\n")
    pms <- normalize.quantiles.use.target(pms, input.vecs$normVec)
  }

  if(verbose) message("Summarizing ...\n")
  pms <- log2(pms)
  
  if(summarize == "average"){  
    exprs <- subColSummarizeAvg((pms-input.vecs$probeVec), pns)
    weights <- NULL
    stderr <- NULL
  }

  if(summarize == "median"){  
    exprs <- subColSummarizeMedian((pms-input.vecs$probeVec), pns)
    weights <- NULL
    stderr <- NULL
  }
  
  if(summarize == "weighted_average"){  
    w <- 1/(input.vecs$probeVarWithin + input.vecs$probeVarBetween)
    if(any(input.vecs$probeVarWithin==0) | any(input.vecs$probeVarBetween==0)) message("Either probeVarWithin or probeVarBetween is 0 for some probes -- setting corresponding weights to 1")
    w[w==Inf] <- 1
    exprs <- subColSummarizeAvg((pms-input.vecs$probeVec)*w, pns)
    W <- as.vector(rowsum(w, pns, reorder=FALSE))
    exprs <- (exprs/W)*as.vector(rowsum(rep(1,length(pns)),pns,reorder=FALSE))
    weights <- NULL
    stderr <- NULL
  }

  if(summarize == "robust_weighted_average"){
    w <- 1/(input.vecs$probeVarWithin + input.vecs$probeVarBetween)
    if(any(input.vecs$probeVarWithin==0) | any(input.vecs$probeVarBetween==0)) message("Either probeVarWithin or probeVarBetween is 0 for some probes -- setting corresponding weights to 1")
    w[w==Inf] <- 1
    N <- 1:dim(pms)[1]
    S <- split(N, pns)
    fit <- lapply(1:length(S), function(i) {
	s <- S[[i]]
	rwaFit2(pms[s,, drop=FALSE], w[s], input.vecs$probeVec[s], input.vecs$probesetSD[i]) 
    })
    names(fit) <- unique(pns)
    exprs <- matrix(unlist(lapply(fit, function(x) x$Estimates)), ncol=ncol(pms), byrow=TRUE)
    rownames(exprs) <- names(fit)
    colnames(exprs) <- colnames(pms)
    if("weights" %in% output.param){
      weights <- matrix(unlist(lapply(fit, function(x) x$Weights)), ncol=ncol(pms), byrow=TRUE)
      rownames(weights) <- pns
      colnames(weights) <- sampleNames(object)
    } else weights <- NULL
    stderr <- matrix(unlist(lapply(fit, function(x) x$StdErrors)), ncol=ncol(pms), byrow=TRUE)
    rownames(stderr) <- names(fit)
    colnames(stderr) <- sampleNames(object)
  }

  if(class(object)!="AffyBatch" & target!="probeset"){
    if(target=="core"){
      featureInfo <- getFidMetaProbesetCore(object)
      exonvec <- input.vecs$probeVecCore
      map <- input.vecs$mapCore
    }
    if(target=="full"){
      featureInfo <- getFidMetaProbesetFull(object)
      exonvec <- input.vecs$probeVecFull
      map <- input.vecs$mapFull
    }
    if(target=="extended"){
      featureInfo <- getFidMetaProbesetExtended(object)
      exonvec <- input.vecs$probeVecExt
      map <- input.vecs$mapExt
    }
    exprs.tmp <- exprs[match(map[,2],rownames(exprs)),,drop=FALSE]
    pns <- map[,1]

    if(summarize == "average" | summarize == "weighted_average"){  
      exprs <- subColSummarizeAvg(exprs.tmp - exonvec, pns)
      weights <- NULL
      stderr <- NULL
    }

    if(summarize == "median"){  
      exprs <- subColSummarizeMedian(exprs.tmp - exonvec, pns)
      weights <- NULL
      stderr <- NULL
    }

    if(summarize == "robust_weighted_average"){
      N <- 1:dim(exprs.tmp)[1]
      S <- split(N, pns)
      fit <- lapply(1:length(S), function(i) {
	s <- S[[i]]
        evec.tmp <- exonvec[s]
        evec.tmp[1] <- evec.tmp[1]-sum(evec.tmp)
	rcModelPLM(exprs.tmp[s,, drop=FALSE], evec.tmp)
      })
      names(fit) <- unique(pns)
      exprs <- matrix(unlist(lapply(fit, function(x) x$Estimates)), ncol=ncol(exprs.tmp), byrow=TRUE)
      rownames(exprs) <- names(fit)
      colnames(exprs) <- colnames(exprs.tmp)
      if("weights" %in% output.param){
        weights <- matrix(unlist(lapply(fit, function(x) x$Weights)), ncol=ncol(exprs.tmp), byrow=TRUE)
        rownames(weights) <- pns
        colnames(weights) <- sampleNames(object)
      } else weights <- NULL
      stderr <- matrix(unlist(lapply(fit, function(x) x$StdErrors)), ncol=ncol(exprs.tmp), byrow=TRUE)
      rownames(stderr) <- names(fit)
      colnames(stderr) <- sampleNames(object)
    }
    if("residuals" %in% output.param){
      residuals <- apply(exprs.tmp,2,function(x) rep(x, table(pns)))
      residuals <- (exprs.tmp - exonvec) - residuals
      rownames(residuals) <- pns
      colnames(residuals) <- sampleNames(object)
    } else residuals <- NULL
  } else {
    if("residuals" %in% output.param){
      residuals <- apply(exprs,2,function(x) rep(x, table(pns)))
      residuals <- (pms-input.vecs$probeVec) - residuals
      rownames(residuals) <- pns
      colnames(residuals) <- sampleNames(object)
    } else residuals <- NULL
  }
  
  colnames(exprs) <- sampleNames(object)
  
  return(list(exprs=exprs, stderr=stderr, weights=weights, residuals=residuals, gammas=NULL))
}

rwaFit2 <- function(x1, x2, x3, x4=NULL){
  ncols <- ncol(x1)
  w.tmp <- x2/max(x2)
  w.tmp <- matrix(rep(w.tmp, ncols), ncol=ncols)
  pe.tmp <- x3
  pe.tmp[1] <- pe.tmp[1]-sum(pe.tmp)
  rcModelWPLM(y=x1, w=w.tmp, row.effects=pe.tmp, input.scale=x4)
}

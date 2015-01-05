frmaFeatureSet <- function(object, background, normalize, summarize, target, input.vecs, output.param, verbose){

  if(background == "rma"){
    if(verbose) message("Background Correcting ...\n")
    object <- backgroundCorrect(object, verbose=verbose)
  }
  
  featureInfo <- getFidProbeset(object)
  pmi <- featureInfo[["fid"]]
  pns <- as.character(featureInfo[["fsetid"]])
  pms <- exprs(object)[pmi,, drop=FALSE]
  
  if(normalize == "quantile"){
    if(verbose) message("Normalizing ...\n")
    pms <- normalize.quantiles.use.target(pms, input.vecs$normVec)
  }
  
  ## first summarize at the probeset-level
  if(verbose) message("Summarizing...\n")
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
    exprs <- matrix(unlist(lapply(fit, function(x) x$Estimates)), ncol=ncol(pms), byrow=TRUE)
    rownames(exprs) <- names(S)
    colnames(exprs) <- sampleNames(object)
    if("weights" %in% output.param){
      weights <- matrix(unlist(lapply(fit, function(x) x$Weights)), ncol=ncol(pms), byrow=TRUE)
      rownames(weights) <- pmi
      colnames(weights) <- sampleNames(object)
    } else weights <- NULL
    stderr <- matrix(unlist(lapply(fit, function(x) x$StdErrors)), ncol=ncol(pms), byrow=TRUE)
    rownames(stderr) <- names(S)
    colnames(stderr) <- sampleNames(object)
  }
  
  ## if selected, use the probeset-level summaries to summarize to another level
  if(target!="probeset"){    
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
      residuals <- apply(exprs,2,function(x) rep(x, table(pns)))
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
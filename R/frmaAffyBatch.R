frmaAffyBatch <- function(object, background, normalize, summarize, input.vecs, output.param, verbose){
  
  if(background == "rma"){
    if(verbose) message("Background Correcting ...\n")
    object <- bg.correct.rma(object)
  }
  
  pms <- pm(object)
  pns <- probeNames(object)
  pmi <- unlist(pmindex(object))
  
  if(!identical(as.character(pmi), names(input.vecs$normVec))){
    map <- match(as.character(pmi), names(input.vecs$normVec))
    input.vecs$normVec <- input.vecs$normVec[map]
    input.vecs$probeVec <- input.vecs$probeVec[map]
    input.vecs$probeVarWithin <- input.vecs$probeVarWithin[map]
    input.vecs$probeVarBetween <- input.vecs$probeVarBetween[map]
    if(!identical(as.character(pmi), names(input.vecs$normVec))){
      stop("Mismatch between pmindex(object) and names of input.vecs and unable to create unique mapping.")
    }
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
    gammas <- NULL
  }
  
  if(summarize %in% c("median","median_polish")){
    if(summarize == "median_polish") message("Note: median and median_polish summarization methods are identical in fRMA.\n")
    exprs <- subColSummarizeMedian((pms-input.vecs$probeVec), pns)
    weights <- NULL
    stderr <- NULL
    gammas <- NULL
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
    gammas <- NULL
  }
  
  if(summarize == "robust_weighted_average"){
    w <- 1/(input.vecs$probeVarWithin + input.vecs$probeVarBetween)
    if(any(input.vecs$probeVarWithin==0) | any(input.vecs$probeVarBetween==0)) message("Either probeVarWithin or probeVarBetween is 0 for some probes -- setting corresponding weights to 1")
    w[w==Inf] <- 1
    N <- 1:dim(pms)[1]
    S <- split(N, pns)
    fit <- lapply(1:length(S), function(i) {
      s <- S[[i]]
      rwaFit2(pms[s,, drop=FALSE], w[s], input.vecs$probeVec[s], input.vecs$probesetSD[names(S)[i]]) 
    })
    
    stderr <- matrix(unlist(lapply(fit, function(x) x$StdErrors)), ncol=ncol(pms), byrow=TRUE)
    rownames(stderr) <- names(S)
    colnames(stderr) <- sampleNames(object)
    
    gammas <- NULL
  }
  
  if(summarize == "random_effect"){
    N <- 1:dim(pms)[1]
    S <- split(N, pns)
    fit <- lapply(1:length(S), function(i) {
      s <- S[[i]]
      batchFit(pms[s,], input.vecs$probeVec[s], input.vecs$probeVarWithin[s], input.vecs$probeVarBetween[s])
    })
    
    if("random_effects" %in% output.param){
      gammas <- matrix(unlist(lapply(fit, function(x) x$Gamma)), ncol=1)
      rownames(gammas) <- pmi
      colnames(gammas) <- NULL
    } else gammas <- NULL
    
    stderr <- NULL
  }
    
  if(summarize %in% c("robust_weighted_average","random_effect")){
    exprs <- matrix(unlist(lapply(fit, function(x) x$Estimates)), ncol=ncol(pms), byrow=TRUE)
    rownames(exprs) <- names(S)
    colnames(exprs) <- sampleNames(object)
    
    if("weights" %in% output.param){
      weights <- matrix(unlist(lapply(fit, function(x) x$Weights)), ncol=ncol(pms), byrow=TRUE)
      rownames(weights) <- pmi
      colnames(weights) <- sampleNames(object)
    } else weights <- NULL
  }
  
  if("residuals" %in% output.param){
    residuals <- apply(exprs, 2, function(x) rep(x, table(pns)))
    residuals <- (pms-input.vecs$probeVec) - residuals
    rownames(residuals) <- pmi
    colnames(residuals) <- sampleNames(object)
  } else residuals <- NULL
  
  return(list(exprs=exprs, stderr=stderr, weights=weights, residuals=residuals, gammas=gammas))
}

frmaRobReg <- function(object, background, normalize, summarize, input.vecs, output.param, verbose){

  cdfname <- cleancdfname(cdfName(object))
  platform <- gsub("cdf","",cdfname)
  
  if(background == "rma"){
    if(verbose) message("Background Correcting ...\n")
    object <- bg.correct.rma(object)
  }
  
  if(is.null(input.vecs$normVec) | is.null(input.vecs$probeVec) | is.null(input.vecs$probeVarWithin) | is.null(input.vecs$probeVarBetween) | (summarize=="robust_weighted_average" & is.null(input.vecs$probesetSD))){
    pkg <- paste(platform, "frmavecs", sep="")
    require(pkg, character.only=TRUE, quiet=TRUE) || stop(paste(pkg, "package must be installed first"))
    data(list=eval(pkg))

    if(is.null(input.vecs$normVec)) input.vecs$normVec <- get(pkg)$normVec
    if(is.null(input.vecs$probeVec)) input.vecs$probeVec <- get(pkg)$probeVec
    if(is.null(input.vecs$probeVarWithin)) input.vecs$probeVarWithin <- get(pkg)$probeVarWithin
    if(is.null(input.vecs$probeVarBetween)) input.vecs$probeVarBetween <- get(pkg)$probeVarBetween
    if(is.null(input.vecs$probesetSD) & summarize=="robust_weighted_average") input.vecs$probesetSD <- get(pkg)$probesetSD
  }

  if(normalize == "quantile"){
    if(verbose) message("Normalizing ...\n")
    pm(object) <- normalize.quantiles.use.target(pm(object), input.vecs$normVec)
  }

  if(verbose) message("Summarizing ...\n")
  pns <- probeNames(object)
  pms <- log2(pm(object))
  
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
    if(any(input.vecs$probeVarWithin==0) | any(input.vecs$probeVarBetween==0)) warning("Either probeVarWithin or probeVarBetween is 0 for some probes -- setting corresponding weights to 1")
    w[w==Inf] <- 1
    exprs <- subColSummarizeAvg((pms-input.vecs$probeVec)*w, pns)
    W <- as.vector(rowsum(w, pns, reorder=FALSE))
    exprs <- (exprs/W)*as.vector(rowsum(rep(1,length(pns)),pns,reorder=FALSE))
    weights <- NULL
    stderr <- NULL
  }

  if(summarize == "robust_weighted_average"){
    w <- 1/(input.vecs$probeVarWithin + input.vecs$probeVarBetween)
    if(any(input.vecs$probeVarWithin==0) | any(input.vecs$probeVarBetween==0)) warning("Either probeVarWithin or probeVarBetween is 0 for some probes -- setting corresponding weights to 1")
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
    if("stderr" %in% output.param){
      stderr <- matrix(unlist(lapply(fit, function(x) x$StdErrors)), ncol=ncol(pms), byrow=TRUE)
      rownames(stderr) <- names(fit)
      colnames(stderr) <- sampleNames(object)
    } else stderr <- NULL
  }

  if("residuals" %in% output.param){
    residuals <- apply(exprs,2,function(x) rep(x, table(pns)))
    residuals <- (pms-input.vecs$probeVec) - residuals
    rownames(residuals) <- pns
    colnames(residuals) <- sampleNames(object)
  } else residuals <- NULL
  
  colnames(exprs) <- sampleNames(object)
  
  return(list(exprs=exprs, stderr=stderr, weights=weights, residuals=residuals))
}

rwaFit2 <- function(x1, x2, x3, x4){
  ncols <- ncol(x1)
  w.tmp <- x2/max(x2)
  w.tmp <- matrix(rep(w.tmp, ncols), ncol=ncols)
  pe.tmp <- x3
  pe.tmp[1] <- pe.tmp[1]-sum(pe.tmp)
  rcModelWPLM(y=x1, w=w.tmp, row.effects=pe.tmp, input.scale=x4)
}


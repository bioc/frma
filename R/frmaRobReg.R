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
    exprs <- subColSummarizeAvg((pms-input.vecs$probeVec)*w, pns)
    W <- as.vector(rowsum(w, pns, reorder=FALSE))
    exprs <- (exprs/W)*as.vector(rowsum(rep(1,length(pns)),pns,reorder=FALSE))
    weights <- NULL
    stderr <- NULL
  }

  if(summarize == "robust_weighted_average"){
    w <- 1/(input.vecs$probeVarWithin + input.vecs$probeVarBetween)
    tmp <- split(data.frame(pms, w, input.vecs$probeVec, rep(input.vecs$probesetSD, table(pns))), pns)
    fit <- lapply(tmp, rwaFit)
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

rwaFit <- function(x){
  pms.tmp <- as.matrix(x[,1:(ncol(x)-3)])
  w.tmp <- x$w/max(x$w)
  w.tmp <- matrix(rep(w.tmp,ncol(pms.tmp)), ncol=ncol(pms.tmp))
  pe.tmp <- x$input.vecs.probeVec
  pe.tmp[1] <- pe.tmp[1]-sum(pe.tmp)
  psd.tmp <- x$probesetSD[1]
  rcModelWPLM(y=pms.tmp, w=w.tmp, row.effects=pe.tmp, input.scale=psd.tmp)
}

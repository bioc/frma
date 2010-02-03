frmaMedPol <- function(object, background, normalize, input.vecs, verbose){

  cdfname <- cleancdfname(cdfName(object))
  platform <- gsub("cdf","",cdfname)

  if(background == "rma"){
    if(verbose) message("Background Correcting ...\n")
    object <- bg.correct.rma(object)
  }
  
  if(is.null(input.vecs$normVec) | is.null(input.vecs$probeVec)){
    pkg <- paste(platform, "frmavecs", sep="")
    require(pkg, character.only=TRUE, quiet=TRUE) || stop(paste(pkg, "package must be installed first"))
    data(list=eval(pkg))

    if(is.null(input.vecs$normVec)) input.vecs$normVec <- get(pkg)$normVec
    if(is.null(input.vecs$probeVec)) input.vecs$probeVec <- get(pkg)$probeVec
  }
  dat <- list(pms=pm(object), pns=probeNames(object), normVec=input.vecs$normVec, probeVec=input.vecs$probeVec)
  
  if(normalize == "quantile"){
      if(verbose) message("Normalizing ...\n")
      dat$pms <- normalize.quantiles.use.target(dat$pms, dat$normVec)
    }
  
  if(verbose) message("Summarizing ...\n")
  dat$pms <- log2(dat$pms) - dat$probeVec
  exprs <- subColSummarizeMedian(dat$pms, dat$pns)

  colnames(exprs) <- sampleNames(object)
  return(list(exprs=exprs, stderr=NULL, weights=NULL, residuals=NULL))
}


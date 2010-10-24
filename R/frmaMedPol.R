frmaMedPol <- function(object, background, normalize, target, input.vecs, verbose){

  if(class(object)=="AffyBatch") cdfname <- cleancdfname(cdfName(object))
  if(class(object)=="ExonFeatureSet") cdfname <- annotation(object)
  platform <- gsub("cdf","",cdfname)

  if(background == "rma"){
    if(verbose) message("Background Correcting ...\n")
    if(class(object)=="AffyBatch") object <- bg.correct.rma(object)
    if(class(object)=="ExonFeatureSet") object <- backgroundCorrect(object, verbose=FALSE)
  }

  if(class(object)=="ExonFeatureSet"){
    if(target=="probeset"){
      featureInfo <- getFidProbeset(object)
    }
    if(target=="core"){
      featureInfo <- getFidMetaProbesetCore(object)
    }
    if(target=="full"){
      featureInfo <- getFidMetaProbesetFull(object)
    }
    if(target=="extended"){
      featureInfo <- getFidMetaProbesetExtended(object)
    }

    pmi <- featureInfo[["fid"]]
    pns <- as.character(featureInfo[["fsetid"]])
    pms <- exprs(object)[pmi,, drop=FALSE]
  }
  if(class(object)=="AffyBatch"){
    pms <- pm(object)
    pns <- probeNames(object)
  }
  
  if(is.null(input.vecs$normVec) | is.null(input.vecs$probeVec)){
    pkg <- paste(platform, "frmavecs", sep="")
    require(pkg, character.only=TRUE, quiet=TRUE) || stop(paste(pkg, "package must be installed first"))
    data(list=eval(pkg))

    if(is.null(input.vecs$normVec)) input.vecs$normVec <- get(pkg)$normVec
    if(is.null(input.vecs$probeVec)) input.vecs$probeVec <- get(pkg)$probeVec
  }
  
  if(normalize == "quantile"){
      if(verbose) message("Normalizing ...\n")
      pms <- normalize.quantiles.use.target(pms, input.vecs$normVec)
    }
  
  if(verbose) message("Summarizing ...\n")
  pms <- log2(pms) - input.vecs$probeVec
  exprs <- subColSummarizeMedian(pms, pns)

  colnames(exprs) <- sampleNames(object)
  return(list(exprs=exprs, stderr=NULL, weights=NULL, residuals=NULL))
}


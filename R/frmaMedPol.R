frmaMedPol <- function(object, background, normalize, target, input.vecs, verbose){

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
     ((is.null(input.vecs$probeVecCore)|is.null(input.vecs$mapCore)) & target=="core") |
     ((is.null(input.vecs$probeVecExt)|is.null(input.vecs$mapExt)) & target=="extended") |
     ((is.null(input.vecs$probeVecFull)|is.null(input.vecs$mapFull)) & target=="full") ){
    pkg <- paste(platform, "frmavecs", sep="")
    require(pkg, character.only=TRUE, quietly=TRUE) || stop(paste(pkg, "package must be installed first"))
    data(list=eval(vecdataname))

    if(is.null(input.vecs$normVec)) input.vecs$normVec <- get(vecdataname)$normVec
    if(is.null(input.vecs$probeVec)) input.vecs$probeVec <- get(vecdataname)$probeVec
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
  }
  
  if(normalize == "quantile"){
      if(verbose) message("Normalizing ...\n")
      pms <- normalize.quantiles.use.target(pms, input.vecs$normVec)
    }

  if(verbose) message("Summarizing ...\n")
  pms <- log2(pms) - input.vecs$probeVec
  exprs <- subColSummarizeMedian(pms, pns)
  
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
    exprs.tmp <- exprs[match(map[,2],rownames(exprs)),,drop=FALSE] - exonvec
    pns <- map[,1]
    exprs <- subColSummarizeMedian(exprs.tmp, pns)
  }

  colnames(exprs) <- sampleNames(object)
  rownames(exprs) <- unique(pns)
  
  return(list(exprs=exprs, stderr=NULL, weights=NULL, residuals=NULL, gammas=NULL))
}


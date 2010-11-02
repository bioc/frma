frma <- function(object, background="rma", normalize="quantile", summarize="robust_weighted_average", target="core", input.vecs=list(normVec=NULL, probeVec=NULL, probeVarBetween=NULL, probeVarWithin=NULL, probesetSD=NULL), output.param=NULL, verbose=FALSE){

  if(!class(object) %in% c("AffyBatch", "ExonFeatureSet")) stop("object must be of class AffyBatch or ExonFeatureSet.")

  if(class(object)=="ExonFeatureSet") target <- match.arg(target, c("core", "full", "extended", "probeset"))
  if(verbose & class(object)=="ExonFeatureSet") message(paste("Exon summarization at ", target, " level.\n", sep=""))

  if(!background %in% c("none", "rma")) stop("background must be either none or rma")
  if(!normalize %in% c("none", "quantile")) stop("normalize must be either none or quantile")
  if(!summarize %in% c("median_polish", "average", "median", "weighted_average", "robust_weighted_average", "random_effect")) stop("summarize must be one of: median_polish, average, median, weighted_average, robust_weighted_average, random_effect")

  if(summarize=="random_effect" & (any(input.vecs$probeVarBetween==0) | any(input.vecs$probeVarWithin==0))) stop("If summarize method is random_effect then probeVarBetween and probeVarWithin must be greater than zero for all probes.")

  if(summarize=="random_effect" & class(object)=="ExonFeatureSet") stop("Summarization method random_effect is not implemented for ExonFeatureSet objects.")
  
  if(class(object)=="AffyBatch") cdfname <- cleancdfname(cdfName(object))
  if(class(object)=="ExonFeatureSet") cdfname <- annotation(object)
  platform <- gsub("cdf","",cdfname)

  if(summarize == "median_polish") output <- frmaMedPol(object, background, normalize, target, input.vecs, verbose)
  if(summarize %in% c("average", "median", "weighted_average", "robust_weighted_average")) output <- frmaRobReg(object, background, normalize, summarize, target, input.vecs, output.param, verbose)
  if(summarize == "random_effect") output <- frmaBatch(object, background, normalize, input.vecs, output.param, verbose)

  if(is.null(output.param)){
    if(is.null(output$stderr)){
      e <- new("ExpressionSet", assayData=assayDataNew(exprs=output$exprs), annotation=platform)
    } else{
      e <- new("ExpressionSet", assayData=assayDataNew(exprs=output$exprs, se.exprs=output$stderr), annotation=platform)
    }
  } else {
    if(is.null(output$stderr)){
      e <- new("frmaExpressionSet", assayData=assayDataNew(exprs=output$exprs), annotation=platform)
    } else{
      e <- new("frmaExpressionSet", assayData=assayDataNew(exprs=output$exprs, se.exprs=output$stderr), annotation=platform)
    }
    if("weights" %in% output.param) w <- output$weights else w <- matrix(nrow=0, ncol=0)
    if("residuals" %in% output.param) r <- output$residuals else r <- matrix(nrow=0, ncol=0)
    e@weights <- w
    e@residuals <- r
  }

  return(e)
}

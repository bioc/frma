frma <- function(object, background="rma", normalize="quantile", summarize="robust_weighted_average", input.vecs=list(normVec=NULL, probeVec=NULL, probeVarBetween=NULL, probeVarWithin=NULL, probesetSD=NULL), output.param=NULL, verbose=FALSE){

  if(!is(object, "AffyBatch")) stop(paste("argument is", class(object), "frma requires AffyBatch"))

  if(!background %in% c("none", "rma")) stop("background must be either none or rma")
  if(!normalize %in% c("none", "quantile")) stop("normalize must be either none or quantile")
  if(!summarize %in% c("median_polish", "average", "median", "weighted_average", "robust_weighted_average", "batch")) stop("summarize must be one of: median_polish, average, median, weighted_average, robust_weighted_average, batch")

  cdfname <- cleancdfname(cdfName(object))
  platform <- gsub("cdf","",cdfname)

  if(summarize == "median_polish") output <- frmaMedPol(object, background, normalize, input.vecs, verbose)
  if(summarize %in% c("average", "median", "weighted_average", "robust_weighted_average")) output <- frmaRobReg(object, background, normalize, summarize, input.vecs, output.param, verbose)
  if(summarize == "batch") output <- frmaBatch(object, background, normalize, input.vecs, output.param, verbose)

  if(is.null(output.param)){
    e <- new("ExpressionSet", exprs=output$exprs, annotation=annotation(object))
  } else if("stderr" %in% output.param & length(output.param)==1){
    e <- new("ExpressionSet", assayData=assayDataNew(exprs=output$exprs, se.exprs=output$stderr), annotation=annotation(object))
  } else{
    if("stderr" %in% output.param){
      e <- new("frmaExpressionSet", assayData=assayDataNew(exprs=output$exprs, se.exprs=output$stderr), annotation=annotation(object))
    } else{
      e <- new("frmaExpressionSet", exprs=output$exprs, annotation=annotation(object))
    }
    if("weights" %in% output.param) e@weights <- output$weights
    if("residuals" %in% output.param) e@residuals <- output$residuals
  }
  
  return(e)
}
  

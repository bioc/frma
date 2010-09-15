frma <- function(object, background="rma", normalize="quantile", summarize="robust_weighted_average", input.vecs=list(normVec=NULL, probeVec=NULL, probeVarBetween=NULL, probeVarWithin=NULL, probesetSD=NULL), output.param=NULL, verbose=FALSE){

  if(!is(object, "AffyBatch")) stop(paste("argument is", class(object), "frma requires AffyBatch"))

  if(!background %in% c("none", "rma")) stop("background must be either none or rma")
  if(!normalize %in% c("none", "quantile")) stop("normalize must be either none or quantile")
  if(!summarize %in% c("median_polish", "average", "median", "weighted_average", "robust_weighted_average", "batch")) stop("summarize must be one of: median_polish, average, median, weighted_average, robust_weighted_average, batch")

  if(summarize=="batch" & (any(input.vecs$probeVarBetween==0) | any(input.vecs$probeVarWithin==0))) stop("If summarize method is batch then probeVarBetween and probeVarWithin must be greater than zero for all probes.")
  
  cdfname <- cleancdfname(cdfName(object))
  platform <- gsub("cdf","",cdfname)

  if(summarize == "median_polish") output <- frmaMedPol(object, background, normalize, input.vecs, verbose)
  if(summarize %in% c("average", "median", "weighted_average", "robust_weighted_average")) output <- frmaRobReg(object, background, normalize, summarize, input.vecs, output.param, verbose)
  if(summarize == "batch") output <- frmaBatch(object, background, normalize, input.vecs, output.param, verbose)

  if("stderr" %in% output.param) stderr <- output$stderr else stderr <- matrix(nrow=0, ncol=0)
  if("weights" %in% output.param) w <- output$weights else w <- matrix(nrow=0, ncol=0)
  if("residuals" %in% output.param) r <- output$residuals else r <- matrix(nrow=0, ncol=0)

  R.model <- PLM.designmatrix3(object)
  
  new("PLMset",
      chip.coefs=output$exprs,
      weights=list("PM.weights"=w, "MM.weights"=matrix(nrow=0, ncol=0)),
      se.chip.coefs=stderr,
      residuals=list("PM.resid"=r, "MM.resid"=matrix(nrow=0, ncol=0)),
      cdfName = object@cdfName,
      phenoData = phenoData(object),
      annotation = annotation(object),
      experimentData = experimentData(object),
      nrow= object@nrow,
      ncol= object@ncol,
      narrays=length(object),
      model.description = list("R.model"=R.model))
}

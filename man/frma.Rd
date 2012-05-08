\name{frma}
\alias{frma}
\title{Frozen Robust Multi-Array Analysis}
\description{
  This function preprocesses an AffyBatch, ExonFeatureSet, or GeneFeatureSet object using the fRMA method.
}
\usage{
frma(object, background="rma", normalize="quantile",
    summarize="robust_weighted_average", target="probeset",
    input.vecs=list(normVec=NULL, probeVec=NULL, probeVarBetween=NULL,
    probeVarWithin=NULL, probesetSD=NULL), output.param=NULL,
    verbose=FALSE) 
}
\arguments{
  \item{object}{an AffyBatch, ExonFeatureSet, or GeneFeatureSet}
  \item{background}{type of background correction to perform: either
    "none" or "rma".}
  \item{normalize}{type of normalization to perform: either "none" or
    "quantile".}  
  \item{summarize}{type of summarization to perform: one of
    "median\_polish", "average", "median", "weighted\_average",
    "robust\_weighted\_average", "random\_effect".}
  \item{target}{summarization level for exon and gene arrays. Must be one of:
    probeset, core, full (exon only), extended (exon only).}
  \item{input.vecs}{a list of vectors to be used in preprocessing. If
    necessary elements are NULL, the correct package with the default
    vectors is loaded if it has been installed. These packages are of
    the form: <platform>frmavecs.}
  \item{output.param}{a vector of output elements to return. By default
    only the expression values and standard errors (if applicable) are
    returned. Additional options are: "residuals", "weights", "random_effects"}
  \item{verbose}{logical value. If TRUE then some messages are displayed
    while the function runs.}
}
\value{
  The function returns an ExpressionSet if output.param=NULL or an
  frmaExpressionSet otherwise.
}
\author{Matthew N. McCall}
\examples{
  library(frmaExampleData)
  data(AffyBatchExample)
  object <- frma(AffyBatchExample)
}
\keyword{manip}


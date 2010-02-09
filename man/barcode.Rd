\name{barcode}
\alias{barcode}
\title{Gene Expression Barcode}
\description{
  This function converts expression values produced via fRMA to a gene
  expression barcode.
}
\usage{
  barcode(object, platform=NULL, mu=NULL, tau=NULL, p=0.5, output="binary")
}
\arguments{
  \item{object}{a vector or matrix of expression values or an
    ExpressionSet or frmaExpressionSet produced by frma}
  \item{platform}{the platform of the input data. One of GPL96, GPL570,
    GPL1261. Required if object is a vector or matrix.}
  \item{mu}{the mean of the unexpressed distribution. If NULL then
    precomputed values are used.}
  \item{tau}{the standard deviation of the unexpressed distribution. If
    NULL then precomputed values are used.}
  \item{p}{prior probability that each gene is expressed.}
  \item{output}{the desired values to be returned. Options are: p-value,
    z-score, binary, or weight.}
}
\value{
  A matrix containing the type of output specified by the output
  parameter. The option \emph{binary} creates a gene expression barcode
  where 1s denote expressed genes and 0s denote unexpressed genes. The
  option \emph{p-value} returns the p-values for the expression values
  under the unexpressed distribution. The option \emph{z-score} returns
  the z-scores for the expression values under the unexpressed
  distribution. The option \emph{weight} returns weights which roughly
  correspond to the probability of expression for each gene. 
}
\author{Matthew N. McCall}
\examples{
  library(frma)
  library(frmaExampleData)
  data(AffyBatchExample)
  object <- frma(AffyBatchExample)
  bc <- barcode(object)
}
\keyword{manip}


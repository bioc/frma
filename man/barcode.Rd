\name{barcode}
\alias{barcode}
\title{Gene Expression Barcode}
\description{
  This function converts expression values produced via fRMA to a gene
  expression barcode.
}
\usage{
  barcode(object, platform=NULL, mu=NULL, tau=NULL, cutoff=6.5, output="binary")
}
\arguments{
  \item{object}{a vector or matrix of expression values or an
    ExpressionSet or frmaExpressionSet produced by frma}
  \item{platform}{the platform of the input data. One of GPL96, GPL570,
    GPL1261. Required if object is a vector or matrix and either mu or
    tau is NULL.}
  \item{mu}{the mean of the unexpressed distribution. If NULL then
    precomputed values are used if possible.}
  \item{tau}{the standard deviation of the unexpressed distribution. If
    NULL then precomputed values are used if possible.}
  \item{cutoff}{the lod score cutoff used if output is binary.}
  \item{output}{the desired values to be returned. Options are: p-value,
    z-score, lod, or binary.}
}
\value{
  A matrix containing the type of output specified by the output
  parameter. The option \emph{binary} creates a gene expression barcode
  where 1s denote expressed genes and 0s denote unexpressed genes. The
  option \emph{p-value} returns the p-values for the expression values
  under the unexpressed distribution. The option \emph{lod} returns the
  LOD scores for expression values under the unexpressed
  distribution. The option \emph{z-score} returns the z-scores for the
  expression values under the unexpressed distribution. 
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


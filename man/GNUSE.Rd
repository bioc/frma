\name{GNUSE}
\alias{GNUSE}
\title{GNUSE}
\description{
  Computes the generalized normalized unscaled standard error (a measure
  of microarray quality).
}
\usage{
  GNUSE(object, medianSE=NULL, type=c("plot", "values", "stats", "density"), ...)
}
\arguments{
  \item{object}{an ExpressionSet or frmaExpressionSet containing
    standard errors produced by frma}
  \item{medianSE}{median standard errors to be used. If NULL, these are
    obtained from the correct frmavecs package.}
  \item{type}{the desired output}
  \item{...}{additional graphical parameters for types plot or density}
}
\value{
  If type is plot, boxplots of GNUSE values are displayed. If type is
  values, the GNUSE values are returned. If type is stats, the median,
  IQR, 95th, and 99th percentiles are reported. If type is density, a
  density plots of GNUSE values are displayed.
}
\author{Matthew N. McCall}
\examples{
  library(frma)
  library(frmaExampleData)
  data(AffyBatchExample)
  object <- frma(AffyBatchExample)
  GNUSE(object, type="stats")
}
\keyword{manip}


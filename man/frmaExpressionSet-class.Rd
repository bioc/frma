\name{frmaExpressionSet-class}
\docType{class}
\alias{class:frmaExpressionSet}
\alias{frmaExpressionSet-class}
\alias{frmaExpressionSet}
\alias{frmaExpressionSet-methods}
\alias{initialize,frmaExpressionSet-method}
\alias{se.exprs,frmaExpressionSet-method}
\alias{se.exprs,ExpressionSet-method}
\alias{weights,frmaExpressionSet-method}
\alias{residuals,frmaExpressionSet-method}
\alias{as.ExpressionSet,frmaExpressionSet-method}
\alias{se.exprs}
\alias{weights}
\alias{residuals}
\alias{as.ExpressionSet}
\title{
  Class to Contain and Describe High-Throughput Expression Level Assays
  preprocessed with fRMA}
\description{This is a class representation for fRMA-preprocessed expression
  data. \code{frmaExpressionSet} class is derived from
  \code{ExpressionSet}, and requires a matrix named \code{exprs}
  and optionally matrices named \code{se.exprs}, \code{weights}, and \code{residuals}.}
\section{Extends}{
  Extends class \code{ExpressionSet}.
  }
\section{Creating Objects}{

  \code{new("frmaExpressionSet",
    phenoData = new("AnnotatedDataFrame"),
    featureData = new("AnnotatedDataFrame"),
    experimentData = new("MIAME"),
    annotation = character(0),
    exprs = new("matrix"),
    weights = new("matrix"),
    residuals = new("matrix")
  }  

  This creates a \code{frmaExpressionSet} with \code{assayData}
  implicitly created to contain \code{exprs}. The only required named
  argument is \code{exprs}. Two optional named matrices, \code{weights}
  and \code{residuals}, can be added to the object.

  \code{new("frmaExpressionSet",
    assayData = assayDataNew(exprs=new("matrix"),se.exprs=new("matrix")),  
    phenoData = new("AnnotatedDataFrame"),
    featureData = new("AnnotatedDataFrame"),
    experimentData = new("MIAME"),
    annotation = character(0),
    weights = new("matrix"),
    residuals = new("matrix")
  }

  This creates a \code{frmaExpressionSet} with \code{assayData}
  provided explicitly. In this form, the only required named argument
  is \code{assayData}. Two optional named matrices, \code{weights}
  and \code{residuals}, can be added to the object.
}
    
\section{Slots}{
  \describe{
    \item{\code{weights}:}{weights used in the summarization step}
    \item{\code{residuals}:}{residuals from fitting the probe-level model}
  }
  Inherited from \code{ExpressionSet}:
  \describe{
    \item{\code{assayData}:}{Contains matrices with equal dimensions,
      and with column number equal to
      \code{nrow(phenoData)}. \code{assayData} must contain a matrix
      \code{exprs} with rows representing features and columns
      representing samples. It may also contain a matrix \code{se.exprs}
    containing standard errors.}  
    \item{\code{phenoData}:}{See \code{eSet}}
    \item{\code{annotation}:}{See \code{eSet}}
    \item{\code{featureData}:}{See \code{eSet}}
    \item{\code{experimentData}:}{See \code{eSet}} 
  }
}

\section{Methods}{
  Class-specific methods:
  \describe{
    \item{\code{se.exprs(frmaExpressionSet)}}{Access elements named
      \code{se.exprs} in the \code{AssayData-class} slot.}
    \item{\code{weights(frmaExpressionSet)}}{Access elements named
      \code{weights}}
    \item{\code{residuals(frmaExpressionSet)}}{Access elements named
      \code{residuals}}
  }

  For derived methods (see \code{ExpressionSet}).
}

\seealso{
  \code{eSet-class}, \code{ExpressionSet-class}, \code{frma}. 
} 

\examples{
# create an instance of frmaExpressionSet
new("frmaExpressionSet")
}

\keyword{classes}
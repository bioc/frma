setMethod("initialize", "frmaExpressionSet",
          function(.Object,
                assayData = assayDataNew(exprs=exprs, se.exprs=se.exprs, ...),
                exprs=new("matrix"),
                se.exprs=new("matrix"),
                weights=new("matrix"),
                residuals=new("matrix"),
                randomeffects=new("matrix"),
                phenoData = annotatedDataFrameFrom  (assayData, byrow=FALSE),
                featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                experimentData = new("MIAME"),
                annotation = new("character"), ...){
            .Object <- callNextMethod(.Object,
                                      assayData = assayData,
                                      phenoData = phenoData,
                                      experimentData = experimentData,
                                      annotation = annotation,
                                      featureData = featureData, ...)
            .Object
          })

## slot getters
setMethod("se.exprs", signature(object="frmaExpressionSet"),
          function(object) assayDataElement(object, "se.exprs"))

weights.frmaExpressionSet <- function(object, ...) object@weights
setMethod("weights", "frmaExpressionSet", weights.frmaExpressionSet)

residuals.frmaExpressionSet <- function(object, ...) object@residuals
setMethod("residuals", "frmaExpressionSet", residuals.frmaExpressionSet)

setMethod("randomeffects", signature(object="frmaExpressionSet"),
          function(object) object@randomeffects)

## validity method
setValidity("frmaExpressionSet", function(object){
  if(!is.null(se.exprs(object))){
    if(any(dim(exprs(object))!=dim(se.exprs(object)))) return("exprs and se.exprs must have the same dimensions")
  }
  if(nrow(weights(object))>0){
    if(ncol(weights(object))!=ncol(exprs(object))) return("weights must have the same number of columns as exprs")
  }
  if(nrow(residuals(object))>0){
    if(ncol(residuals(object))!=ncol(exprs(object))) return("residuals must have the same number of columns as exprs")
  }
  if(nrow(residuals(object))>0 & nrow(weights(object))>0){
    if(nrow(residuals(object)) != nrow(weights(object))) return("if not empty, residuals and weights must have the same number of rows")
  }
  if(nrow(randomeffects(object))>0 & nrow(weights(object))>0){
    if(nrow(randomeffects(object)) != nrow(weights(object))) return("if not empty, randomeffects and weights must have the same number of rows")
  }
  if(nrow(randomeffects(object))>0 & nrow(residuals(object))>0){
    if(nrow(randomeffects(object)) != nrow(residuals(object))) return("if not empty, randomeffects and residuals must have the same number of rows")
  }
  TRUE
})

## slot setters
setReplaceMethod("se.exprs", signature(object="frmaExpressionSet"), function(object, value){
  object <- assayDataElementReplace(object, "se.exprs", value)
  validObject(object)
  object
})

setReplaceMethod("weights", "frmaExpressionSet", function(object, value){
  object@weights <- value
  validObject(object)
  object
})

setReplaceMethod("residuals", "frmaExpressionSet", function(object, value){
  object@residuals <- value
  validObject(object)
  object
})

setReplaceMethod("randomeffects", "frmaExpressionSet", function(object, value){
  object@randomeffects <- value
  validObject(object)
  object
})

## coercion methods
setAs("frmaExpressionSet", "ExpressionSet", function(from){
  new("ExpressionSet", assayData=assayDataNew(exprs=exprs(from), se.exprs=se.exprs(from)), annotation=annotation(from))
})
                                            

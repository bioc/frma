setMethod("initialize", "frmaExpressionSet",
          function(.Object,
                   assayData = assayDataNew(exprs=exprs, ...),
                   exprs=new("matrix"),
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

setMethod("se.exprs", signature(object="frmaExpressionSet"),
          function(object) assayDataElement(object, "se.exprs"))

setMethod("weights", signature(object="frmaExpressionSet"),
          function(object) object@weights)

setMethod("residuals", signature(object="frmaExpressionSet"),
          function(object) object@residuals)

setMethod("as.ExpressionSet", signature(object="frmaExpressionSet"),
          function(object){
            new("ExpressionSet", assayData=assayDataNew(exprs=exprs(object), se.exprs=se.exprs(object)), annotation=annotation(object))
          })
                                            

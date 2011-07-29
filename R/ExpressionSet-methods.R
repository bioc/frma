setMethod("se.exprs", signature(object="ExpressionSet"),
          function(object) assayDataElement(object, "se.exprs"))

setReplaceMethod("se.exprs", signature(object="ExpressionSet"), function(object, value){
  object <- assayDataElementReplace(object, "se.exprs", value)
  validObject(object)
  object
})

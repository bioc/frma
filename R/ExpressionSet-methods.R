setMethod("se.exprs", signature(object="ExpressionSet"),
          function(object) assayDataElement(object, "se.exprs"))

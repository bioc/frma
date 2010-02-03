setClass("frmaExpressionSet", representation(weights="matrix", residuals="matrix"), prototype=list(weights=matrix(), residuals=matrix()), contains="ExpressionSet")


barcode <- function(object, platform=NULL, mu=NULL, tau=NULL, cutoff=6.5, output="binary"){

  if(!class(object) %in% c("matrix", "ExpressionSet", "frmaExpressionSet") & !is.vector(object)) stop("Object must be one of: vector, matrix, ExpressionSet, or frmaExpressionSet.")

  if(!output %in% c("binary", "weight", "z-score", "p-value", "lod")) stop("Output must be one of: binary, weight, z-score, p-value, lod.") 

  if(is.vector(object)) object <- as.matrix(object)
  
  if(class(object) == "frmaExpressionSet"){
    object <- as.ExpressionSet(object)
  }
  
  if(is.matrix(object)){
    if(is.null(platform)) stop("If object is of class matrix, platform cannot be NULL.")
    if(!platform %in% c("GPL96", "GPL570", "GPL1261")) stop("Platform must be one of: GPL96, GPL570, GPL1261")
  } else{
    if(annotation(object) %in% c("hgu133a", "hgu133plus2", "moe430_2")){
      if(cleancdfname(annotation(object)) == "hgu133acdf") platform <- "GPL96"
      if(cleancdfname(annotation(object)) == "hgu133plus2cdf") platform <- "GPL570"
      if(cleancdfname(annotation(object)) == "mouse4302cdf") platform <- "GPL1261"
    } else stop("Microarray platform given by object annotation not recognized -- please supply platform type.")
    object <- as.matrix(exprs(object))
  }

  i.remove <- grep("AFFX", rownames(object))
  if(length(i.remove)>0) object <- as.matrix(object[-i.remove,])

  if(platform == "GPL96" & nrow(object)!=22215) stop("Object does not have the correct dimensions for platform GPL96.")
  if(platform == "GPL570" & nrow(object)!=54613) stop("Object does not have the correct dimensions for platform GPL570.")
  if(platform == "GPL1261" & nrow(object)!=45037) stop("Object does not have the correct dimensions for platform GPL1261.")
  
  if(is.null(mu) | is.null(tau)){
    if(platform=="GPL96") pkg <- "hgu133abarcodevecs"
    if(platform=="GPL570") pkg <- "hgu133plus2barcodevecs"
    if(platform=="GPL1261") pkg <- "mouse4302barcodevecs"
    require(pkg, character.only=TRUE, quiet=TRUE) || stop(paste(pkg, "package must be installed first"))
    data(list=paste("bcparams-", platform, sep=""))
    if(is.null(mu)) mu <- bcparams[,2]
    if(is.null(tau)) tau <- sqrt(bcparams[,3])
  }
  
  e <- object
  num <- round(nrow(e)/10)
  
  if(output %in% c("p-value", "lod", "binary")){
    pval <- pnorm(e, mean=mu, sd=tau, lower.tail=FALSE)
    if(output == "p-value"){
      colnames(pval) <- colnames(e)
      rownames(pval) <- rownames(e)
      return(pval)
    } else{
      lod <- -log10(pval)
      if(output == "lod"){
        colnames(lod) <- colnames(e)
        rownames(lod) <- rownames(e)
        return(lod)
      } else{
        bc <- matrix(as.integer(lod > cutoff), ncol=ncol(e))
        colnames(bc) <- colnames(e)
        rownames(bc) <- rownames(e)
        return(bc)
      }
    }
  }

  if(output == "z-score"){
    z <- sweep(sweep(e, 1, mu), 1, tau, FUN="/")
    colnames(z) <- colnames(e)
    rownames(z) <- rownames(e)
    return(z)
  }
  
  if(output == "weight"){
    unf <- dunif(e, mu, 15)
    nrm <- dnorm(e, mean=mu, sd=tau)
    w <- matrix(ifelse(nrm==0 & unf==0, 0, (p*unf) / ((p*unf) + ((1-p)*nrm))), nrow=nrow(e), ncol=ncol(e))
    colnames(w) <- colnames(e)
    rownames(w) <- rownames(e)
    return(w)
  }
}
